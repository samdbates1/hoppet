!! An example program using a tabulation. It outputs a subset of
!! table 15 of hep-ph/0511119 and this output should be identical
!! to the contents of the file tabulation_example.default_output
!!
!! NB: for the full functionality used in generating the HeraLHC and
!!     Les Houches comparison tables, see ../benchmarks/benchmarks.f90
!!     and carefully read the comments at the beginning. Subtleties
!!     exist in particular wrt the treatment of scales muF/=muR.
!!
!! NB: commented code shows usage with LHAPDF (e.g. for cross-checking 
!!     some public PDF set's evolution) -- to use this part of the code
!!     you will also need to link with LHAPDF
!!
program tabulation_example
     use hoppet
  !! if using LHAPDF, rename a couple of hoppet functions which
  !! would otherwise conflict with LHAPDF 
     !use hoppet, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  implicit none
  real(dp) :: dy, ymax
  integer  :: order, nloop
  !! holds information about the grid
  type(grid_def) :: grid
  !! holds the splitting functions
  type(dglap_holder) :: dh
  !! hold the PDF tabulation
  type(pdf_table)       :: table
  !! hold the coupling
  real(dp)               :: quark_masses(4:6)
  type(running_coupling) :: coupling
  !! hold the initial pdf
  real(dp), pointer :: pdf0(:,:)
  real(dp) :: Q0
  !! hold results at some x, Q
  real(dp) :: Q, pdf_at_xQ(-6:6)
  real(dp), parameter :: heralhc_xvals(9) = &
       & (/1e-5_dp,1e-4_dp,1e-3_dp,1e-2_dp,0.1_dp,0.3_dp,0.5_dp,0.7_dp,0.9_dp/)
  integer  :: ix
  real(dp) :: integral

  !! define the interfaces for LHA pdf (by default not used)
  !! (NB: unfortunately this conflicts with an internal hoppet name,
  !! so make sure that you "redefine" the internal hoppet name, 
  !! as illustrated in the commented "use" line above:
     !! use hoppet, EvolvePDF_hoppet => EvolvePDF, ...)
  ! interface
  !    subroutine EvolvePDF(x,Q,res)
  !      use types; implicit none
  !      real(dp), intent(in)  :: x,Q
  !      real(dp), intent(out) :: res(*)
  !    end subroutine EvolvePDF
  ! end interface

  ! set up parameters for grid
  order = -6
  ymax  = 12.0_dp
  dy    = 0.1_dp

  ! set up the grid itself (this call sets up a nested grid composed of 4 subgrids)
  call InitGridDefDefault(grid, dy, ymax, order=order)

  ! initialise the splitting-function holder
  nloop = 1
  call InitDglapHolder(grid,dh,factscheme=factscheme_FragMSbar,&
       &                      nloop=nloop,nflo=3,nfhi=6)
  write(6,'(a)') "Splitting functions initialised!"

  !! set up LHAPDF with cteq
  !call InitPDFsetByName("cteq61.LHgrid")
  !call InitPDF(0)
  ! allocate and set up the initial pdf from LHAPDF ...
  !Q0 = 5.0_dp 
  !call AllocPDF(grid, pdf0)
  !call InitPDF_LHAPDF(grid, pdf0, EvolvePDF, Q0)

  ! initialise a PDF from the function below (must be contained,
  ! in a "used" module, or with an explicitly defined interface)
  call AllocPDF(grid, pdf0)
  pdf0 = unpolarized_dummy_ff(xValues(grid))
  Q0 = 5.0_dp  ! the initial scale

  ! allocate and initialise the running coupling with a given
  ! set of quark masses (NB: charm mass just above Q0).
  quark_masses(4:6) = (/1.414213563_dp, 4.5_dp, 175.0_dp/)
  call InitRunningCoupling(coupling,alfas=0.35_dp,Q=Q0,nloop=nloop,&
       &                   quark_masses = quark_masses)

  ! create the tables that will contain our copy of the user's pdf
  ! as well as the convolutions with the pdf.
  call AllocPdfTable(grid, table, Qmin=1.0_dp, Qmax=10000.0_dp, & 
       & dlnlnQ = dy/4.0_dp, freeze_at_Qmin=.true.)
  ! add information about the nf transitions to the table (improves
  ! interpolation quality)
  call AddNfInfoToPdfTable(table,coupling)

  ! create the tabulation based on the evolution of pdf0 from scale Q0
  call EvolvePdfTable(table, Q0, pdf0, dh, coupling, nloop=nloop)
  ! alternatively "pre-evolve" so that subsequent evolutions are faster
  !call PreEvolvePdfTable(table, Q0, dh, coupling)
  !call EvolvePdfTable(table,pdf0)
  write(6,'(a)') "Evolution done!"

  ! get the value of the tabulation at some point
  Q = 100.0_dp
  write(6,'(a)')
  write(6,'(a,f8.3,a)') "           Evaluating PDFs at Q = ",Q," GeV"
  write(6,'(a5,2a12,a14,a10,a12)') "x",&
       & "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon"
  do ix = 1, size(heralhc_xvals)
     call EvalPdfTable_xQ(table,heralhc_xvals(ix),Q,pdf_at_xQ)
     write(6,'(es7.1,5es12.4)') heralhc_xvals(ix), &
          &  pdf_at_xQ(2)-pdf_at_xQ(-2), &
          &  pdf_at_xQ(1)-pdf_at_xQ(-1), &
          &  2*(pdf_at_xQ(-1)+pdf_at_xQ(-2)), &
          &  (pdf_at_xQ(-4)+pdf_at_xQ(4)), &
          &  pdf_at_xQ(0)
  end do

  call EvolvePDF(dh,pdf0,coupling,Q0,Q)
  ! Do the integration
  integral = zero
  do ix = -6, 6
    write (*,*) ix
    write (*,*) TruncatedMoment(grid,pdf0(:,ix),1.0_dp)
    integral = integral + TruncatedMoment(grid,pdf0(:,ix),1.0_dp)
  end do
  write (*,*) integral
  
  
  ! some cleaning up (not strictly speaking needed, but illustrates
  ! how it's done)
  call Delete(table)
  call Delete(pdf0)
  call Delete(dh)
  call Delete(coupling)
  call Delete(grid)

contains 
  !======================================================================
  !! The dummy PDF suggested by Vogt as the initial condition for the 
  !! unpolarized evolution (as used in hep-ph/0511119).
  function unpolarized_dummy_ff(xvals) result(pdf)
    real(dp), intent(in) :: xvals(:)
    real(dp)             :: pdf(size(xvals),ncompmin:ncompmax)
    real(dp)             :: quark(size(xvals)), gluon(size(xvals))
    !---------------------
    real(dp), parameter :: N_g = 5.0_dp, N_q = 6.0_dp ! must be > 1
  
    pdf = zero
    ! clean method for labelling as PDF as being in the human representation
    ! (not actually needed after setting pdf=0
    call LabelPdfAsHuman(pdf)

    gluon = N_g * (N_g-1) * xvals * (1-xvals)**(N_g-2)
    quark = N_q * (N_q-1) * xvals * (1-xvals)**(N_q-2)
    ! labels iflv_g, etc., come from the hoppet module, inherited
    ! from the main program
    ! remember that these are all x*q(x)
    pdf(:, iflv_g) = gluon
    pdf(:, iflv_d) = quark
    pdf(:,-iflv_d) = quark
    pdf(:, iflv_u) = quark
    pdf(:,-iflv_u) = quark
    pdf(:, iflv_s) = quark
    pdf(:,-iflv_s) = quark
    pdf(:, iflv_c) = quark
    pdf(:,-iflv_c) = quark
    pdf(:, iflv_b) = quark
    pdf(:,-iflv_b) = quark

  end function unpolarized_dummy_ff

end program tabulation_example


