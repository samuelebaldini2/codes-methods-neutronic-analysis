subroutine compute_mie(mr, mi, rad, lam,                  &
                       size_parameter,                    &
                       extinction_efficiency,             &
                       scattering_efficiency,             &
                       absorption_efficiency,             &
                       sphere_geometric_cross_section,    &
                       single_scattering_albedo,          &
                       asymmetry_factor,                  &
                       radiation_pressure_efficiency,     &
                       number_of_moments,                 &
                       legendre_moments)

use iso_fortran_env, only: real64
implicit none

! ============================================================
! Inputs
! ============================================================

real, intent(in) :: mr
real, intent(in) :: mi
real, intent(in) :: rad
real, intent(in) :: lam

! ============================================================
! Outputs
! ============================================================

real, intent(out) :: size_parameter
real, intent(out) :: extinction_efficiency
real, intent(out) :: scattering_efficiency
real, intent(out) :: absorption_efficiency
real, intent(out) :: sphere_geometric_cross_section
real, intent(out) :: single_scattering_albedo
real, intent(out) :: asymmetry_factor
real, intent(out) :: radiation_pressure_efficiency

integer, intent(out) :: number_of_moments

! Transfer buffer to Python.
! Only elements 1:number_of_moments+1 are meaningful.
real, intent(out) :: legendre_moments(2205)


interface

    subroutine MIEV0(xx, crefin, perfct, mimcut, anyang,       &
                     numang, xmu, nmom, ipolzn, momdim, prnt,  &
                     qext, qsca, gqsc, pmom, sforw, sback, s1, &
                     s2, tforw, tback, spike)

        implicit none

        logical :: anyang, perfct, prnt(*)

        integer :: ipolzn
        integer :: momdim
        integer :: numang
        integer :: nmom

        real :: gqsc
        real :: mimcut
        real :: pmom(0:momdim,*)
        real :: qext
        real :: qsca
        real :: spike
        real :: xmu(*)
        real :: xx

        complex :: crefin
        complex :: sforw
        complex :: sback
        complex :: s1(*)
        complex :: s2(*)
        complex :: tforw(*)
        complex :: tback(*)

    end subroutine MIEV0

end interface


! ============================================================
! Constants
! ============================================================

real(kind=real64), parameter :: pi = 3.141592653589793_real64
real(kind=real64), parameter :: two_pi = 2 * pi
real(kind=real64), parameter :: one_third = 1._real64 / 3._real64


! ============================================================
! MIEV0 arguments
! ============================================================

real :: xx
complex :: crefin
logical :: perfct
real :: mimcut
logical :: anyang

integer :: numang

real, allocatable :: xmu(:)

integer :: nmom
integer :: ipolzn
integer :: momdim

logical :: prnt(2)

real :: qext
real :: qsca
real :: gqsc

complex, allocatable :: s1(:)
complex, allocatable :: s2(:)

complex :: sforw
complex :: sback
complex :: tforw(2)
complex :: tback(2)

real :: spike

real, allocatable :: pmom(:, :)


! ============================================================
! Local variables
! ============================================================

integer :: i
integer :: n
integer :: step

real(kind=real64) :: xx_dp
real(kind=real64) :: fnorm


! ============================================================
! Adaptation of mieleg.F90
! ============================================================

crefin = cmplx(mr, mi)

prnt = [.false., .false.]
perfct = .false.
anyang = .false.
mimcut = 1.e-6
ipolzn = 0
numang = 1000001

xx_dp = two_pi * real(rad, kind=real64) / &
         real(lam, kind=real64)

step = 2 / (numang - 1)

nmom = int(2 * (xx_dp + 4 * xx_dp**one_third + 2))

xx = real(xx_dp)

momdim = nmom + 1

allocate(xmu(numang))
allocate(s1(numang))
allocate(s2(numang))
allocate(pmom(0:momdim, 1))

do n = 1, numang
    xmu(n) = real(1 - (n - 1) * step)
end do

call MIEV0(xx, crefin, perfct, mimcut, anyang,       &
           numang, xmu, nmom, ipolzn, momdim, prnt,  &
           qext, qsca, gqsc, pmom, sforw, sback, s1, &
           s2, tforw, tback, spike)

fnorm = 4._real64 / &
        (xx_dp**2 * real(qsca, kind=real64))

size_parameter = xx

extinction_efficiency = qext
scattering_efficiency = qsca

absorption_efficiency = qext - qsca

sphere_geometric_cross_section = &
    real(pi * real(rad, kind=real64)**2)

single_scattering_albedo = qsca / qext

asymmetry_factor = gqsc / qsca

radiation_pressure_efficiency = qext - gqsc

number_of_moments = nmom

legendre_moments = 0.0

do i = 0, nmom

    legendre_moments(i + 1) = &
        real(fnorm * real(pmom(i,1), kind=real64))

end do

deallocate(xmu)
deallocate(s1)
deallocate(s2)
deallocate(pmom)

end subroutine compute_mie
