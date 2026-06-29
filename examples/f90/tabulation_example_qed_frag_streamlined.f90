

!! An example program using HOPPET's streamlined interface to
!! perform a QCD+QED evolution of PDFs.

!! A module that contains the initial condition for the evolution
! module lhasub4streamlinedqed
!   implicit none

! contains
!   !======================================================================
!   !! The dummy PDF suggested by Vogt as the initial condition for the 
!   !! unpolarized evolution (as used in hep-ph/0511119), 
!   !! extended to include QED partons in a simple way.
!   subroutine LHAsub(x,Q,xpdf)
!   use types; use consts_dp; use hoppet; implicit none
!     real(dp), intent(in)  :: x,Q
!     real(dp), intent(out) :: xpdf(-6:ncompmaxLeptons)
!     real(dp) :: gluon, uv, dv
!     real(dp) :: ubar, dbar
!     !---------------------
!     real(dp), parameter :: N_g = 1.7_dp, N_ls = 0.387975_dp
!     real(dp), parameter :: N_uv=5.107200_dp, N_dv = 3.064320_dp
!     real(dp), parameter :: N_db = half*N_ls
  
!     ! Set to zero the xpdf array
!     xpdf = zero
    
  
!     !-- remember that these are all x*q(x)
!     gluon = N_g * x**(-0.1_dp) * (1-x)**5
!     uv = N_uv * x**0.8_dp * (1-x)**3
!     dv = N_dv * x**0.8_dp * (1-x)**4
!     dbar = N_db * x**(-0.1_dp) * (1-x)**6
!     ubar = dbar * (1-x)
  
!     ! labels iflv_g, etc., come from the hoppet module, inherited
!     ! from the main program
!   end subroutine LHAsub
! end module lhasub4streamlinedqed

module lhasub4streamlinedqed
  implicit none

contains
  !======================================================================
  !! dummy fragmentation containing only u and g.
  subroutine LHAsub(x,Q,xpdf)
  use types; use consts_dp; use hoppet; implicit none
    real(dp), intent(in)  :: x,Q
    real(dp), intent(out) :: xpdf(-6:6)
    real(dp) :: u, glu
    real(dp) :: N_u, N_g ! must be > 1
  
    ! Set to zero the xpdf array
    xpdf = zero

    N_u = 5
    N_g = 6
    
    ! labels iflv_g, etc., come from the hoppet module, inherited
    ! from the main program
    !-- remember that these are all x*q(x)
    xpdf(iflv_g) = N_g * (N_g-1) * x * (1-x)**(N_g-2)
    xpdf(iflv_u) = N_u * (N_u-1) * x * (1-x)**(N_u-2)
  end subroutine LHAsub
end module lhasub4streamlinedqed

! ! This module was partly written by GPT-5.5
! module pdf_integrator_module
!   use types
!   use hoppet
!   use integrator
!   implicit none

!   type :: pdf_integrator
!     real(dp) :: xmin = 0.0_dp    ! set with xmin = exp(-ymax)
!     real(dp) :: EPS  = 1.0e-8_dp ! precision with which to integrate.
!   contains
!     procedure :: integrate
!   end type pdf_integrator

!   type, extends(ignd_class) :: pdf_integrand
!     integer  :: iflavour
!     real(dp) :: Q
!     logical  :: divide_by_x = .false.
!   contains
!     procedure :: f => pdf_integrand_f
!   end type pdf_integrand

! contains

!   ! Returns the PDF as a function of x only
!   ! for feeding to the integrator
!   function pdf_integrand_f(this, x) result(val)
!     class(pdf_integrand), intent(in) :: this
!     real(dp), intent(in) :: x
!     real(dp)             :: val
!     real(dp)             :: pdf(-6:ncompmaxLeptons)

!     call hoppetEval(x, this%Q, pdf)

!     val = pdf(this%iflavour)

!     ! Optional: return f(x) instead of x*f(x)
!     if (this%divide_by_x) then
!       if (x > 0.0_dp) then
!         val = val / x
!       else
!         val = 0.0_dp
!       endif 
!     endif

!   end function pdf_integrand_f

!   ! integrates the pdf of a given iflav at a given Q.
!   ! Hoppet pdfs are stored as x*f(x), so to evaluate
!   ! the integral of f(x), 'divide_by_x' must be true.
!   function integrate(this, iflav, Q, divide_by_x) result(integral)
!     class(pdf_integrator), intent(in) :: this
!     integer,  intent(in) :: iflav
!     real(dp), intent(in) :: Q
!     logical,  intent(in) :: divide_by_x
!     real(dp) :: integral
!     type (pdf_integrand) :: pdf

!     pdf%iflavour    = iflav
!     pdf%Q           = Q
!     pdf%divide_by_x = divide_by_x

!     integral = ig_LinWeight( &
!                 pdf, &
!                 this%xmin, 1.0_dp, &
!                 1.0_dp, 1.0_dp, &
!                 this%EPS)
!   end function integrate


! end module pdf_integrator_module


!! the main program illustrating the streamlined interface with QED
program tabulation_example_qed_streamlined
  use hoppet; use lhasub4streamlinedqed; use convolution
  !! if using LHAPDF, rename a couple of hoppet functions which
  !! would otherwise conflict with LHAPDF 
  !use hoppet, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  implicit none
  real(dp) :: dy, ymax, dlnlnQ, Qmin, Qmax, muR_Q
  real(dp) :: xmin
  real(dp) :: asQ, Q0alphas, Q0pdf
  real(dp) :: mc,mb,mt
  integer  :: order, nloop
  !! hold results at some x, Q
  real(dp) :: Q, xpdf_at_xQ(-6:ncompmaxLeptons)
  real(dp), parameter :: heralhc_xvals(9) = &
       & (/1e-5_dp,1e-4_dp,1e-3_dp,1e-2_dp,0.1_dp,0.3_dp,0.5_dp,0.7_dp,0.9_dp/)
  integer  :: ix
  logical  :: use_qed, use_qcd_qed, use_Plq_nnlo
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


  ! Streamlined initialization
  ! including  parameters for x-grid
  order = -6
  ymax  = 12.0_dp
  dy    = 0.025_dp
  xmin  = exp(-ymax)
  ! and parameters for Q tabulation
  Qmin  = 1.0_dp
  Qmax  = 28000.0_dp
  dlnlnQ = dy/4.0_dp
  ! and number of loops to initialise!
  nloop = 1

  ! use_qed = .true.
  ! use_qcd_qed  = .true.
  ! use_Plq_nnlo = .false.
  ! call hoppetSetQED(use_qed, use_qcd_qed, use_Plq_nnlo)
  ! call with_qed_true
  ! call with_qcd_qed_true
  ! call with_Plq_false 
  
  call hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,&
       &         order,factscheme_FragMSbar)
  write(6,'(a)') "Streamlined initialization completed!"
  
  ! Set heavy flavour scheme
  mc = 1.414213563_dp   ! sqrt(2.0_dp) + epsilon
  mb = 4.5_dp
  mt = 175.0_dp
  call hoppetSetVFN(mc, mb, mt)

  ! Streamlined evolution
 
  ! Set parameters of running coupling
  asQ = 0.35_dp
  Q0alphas = sqrt(2.0_dp)
  muR_Q = 1.0_dp
  Q0pdf = 5.0_dp ! The initial evolution scale

  call hoppetEvolve(asQ, Q0alphas, nloop,muR_Q, LHAsub, Q0pdf)

  ! get the value of the tabulation at some point
  Q = 10.0_dp
  write(6,'(a)')
  write(6,'(a,f8.3,a)') "           Evaluating PDFs at Q = ",Q," GeV"
  write(6,'(a5,2a12,a14,a11,5a12)') "x",&
       & "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon",&
       & "photon","e+ + e-","mu+ + mu-"," tau+ + tau-"
  do ix = 1, size(heralhc_xvals)
     call hoppetEval(heralhc_xvals(ix),Q,xpdf_at_xQ)
     write(6,'(es7.1,9es12.4)') heralhc_xvals(ix), &
          &  xpdf_at_xQ(2)-xpdf_at_xQ(-2), &
          &  xpdf_at_xQ(1)-xpdf_at_xQ(-1), &
          &  2*(xpdf_at_xQ(-1)+xpdf_at_xQ(-2)), &
          &  (xpdf_at_xQ(-4)+xpdf_at_xQ(4)), &
          &  xpdf_at_xQ(0),  &
          &  xpdf_at_xQ(iflv_photon), &
          &  xpdf_at_xQ(iflv_electron), &
          &  xpdf_at_xQ(iflv_muon),&
          &  xpdf_at_xQ(iflv_tau)
  end do

  ! NB: there is no cached evolution option available
  ! with QED turned on

  !call hoppetWriteLHAPDFGrid("test_qed",0)
  
  ! perform cleanup (not strictly required)
  call hoppetDeleteAll()

end program tabulation_example_qed_streamlined
