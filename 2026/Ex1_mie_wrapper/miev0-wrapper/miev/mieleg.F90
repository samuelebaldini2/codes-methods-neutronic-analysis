! SPDX-FileCopyrightText: 2026 Alberto P
!
! SPDX-License-Identifier: MPL-2.0

! Compute Mie scattering properties for a single sphere using the MIEV0 code.
program miescat
use iso_fortran_env, only: real64
implicit none(type, external)

! interface to MIEV0 subroutine
interface
    subroutine MIEV0(xx, crefin, perfct, mimcut, anyang,       &
                     numang, xmu, nmom, ipolzn, momdim, prnt,  &
                     qext, qsca, gqsc, pmom, sforw, sback, s1, &
                     s2, tforw, tback, spike)
        implicit none(type)
        logical :: anyang, perfct, prnt(*)
        integer :: ipolzn, momdim, numang, nmom
        real    :: gqsc, mimcut, pmom(0:momdim,*), qext, qsca, spike, &
                   xmu(*), xx
        complex :: crefin, sforw, sback, s1(*), s2(*), tforw(*), tback(*)
    end subroutine MIEV0
end interface

! Global constants
real(kind=real64), parameter :: pi = 3.141592653589793_real64
real(kind=real64), parameter :: two_pi = 2 * pi
real(kind=real64), parameter :: one_third = 1._real64 / 3._real64

! MIEV0 arguments
! inputs
real :: xx                            ! Mie size parameter ( 2 * pi * radius / wavelength )
complex :: crefin                     ! Refractive index of the particle (mr  j*mi)
logical :: perfct                     ! True if refractive index is infinite
real :: mimcut                        ! minimum imaginary refractive index
logical :: anyang                     ! symmetric or not angles
integer :: numang                     ! number of angles for scattering amplitude
real, allocatable :: xmu(:)           ! cosine of angles for scattering amplitude functions
integer :: nmom                       ! highest Legendre moments to compute
integer :: ipolzn                     ! polarization calculation type (0, 1, or 2)
integer :: momdim                     ! first dimension of the array of Legendre moments
logical :: prnt(2)                    ! print flags for S1/S2 (pos 1) and other diagnostic info (pos 2)
! outputs
real :: qext                          ! extinction efficiency factor
real :: qsca                          ! scattering efficiency factor
real :: gqsc                          ! asymmetry parameter times scattering efficiency
complex, allocatable :: s1(:), s2(:)  ! scattering amplitude functions at angles specified in xmu
complex :: sforw                      ! forward scattering amplitude
complex :: sback                      ! backward scattering amplitude
complex :: tforw(2)                   ! forward scattering values for T1 and T2
complex :: tback(2)                   ! backward scattering values for T1 and T2
real :: spike                         ! magnitude of the smallest denominator of Mie coefficients
real, allocatable :: pmom(:, :)       ! Legendre moments of the scattering phase function (1st index is moment order, 2nd index is polarization type)

! command-line arguments
integer :: num_args
character(len=32) :: arg
real :: arg_real

! local indexes, error states, and support variables
integer :: i, n
integer :: errstat
integer :: step
real(kind=real64) :: xx_dp

! input parameters
real :: mr   ! real part of the refractive index
real :: mi   ! imaginary part of the refractive index
real :: rad  ! particle radius
real :: lam  ! wavelength of incident light

! Mie scattering variables
real(kind=real64) :: fnorm    ! Wiscombe normalization factor for Legendre moments


! read command line arguments
num_args = command_argument_count()
if (num_args /= 4) then
    write(*, *) "Usage: miescat <mr> <mi> <rad> <lam>"
    stop 1
end if

call get_command_argument(1, arg)
read(arg, *, iostat=errstat) arg_real
if (errstat /= 0) then
    write(*, *) "Error: Argument is not a valid real number."
    stop 2
end if
mr = real(arg_real)

call get_command_argument(2, arg)
read(arg, *, iostat=errstat) arg_real
if (errstat /= 0) then
    write(*, *) "Error: Argument is not a valid real number."
    stop 2
end if
mi = real(arg_real)

crefin = cmplx(mr, mi)

call get_command_argument(3, arg)
read(arg, *, iostat=errstat) arg_real
if (errstat /= 0) then
    write(*, *) "Error: Argument is not a valid real number."
    stop 2
end if
rad = real(arg_real)

call get_command_argument(4, arg)
read(arg, *, iostat=errstat) arg_real
if (errstat /= 0) then
    write(*, *) "Error: Argument is not a valid real number."
    stop 2
end if
lam = real(arg_real)

! set computation flags for MIEV0 calculation
#ifndef NDEBUG
prnt = [.false., .true.]
#else
prnt = [.false., .false.]
#endif
perfct = .false.
anyang = .false.
mimcut = 1.e-6
ipolzn = 0
numang = 1000001

! size parameter for Mie calculation, computed in double precision to avoid precision loss
xx_dp = two_pi * real(rad, kind=real64) / real(lam, kind=real64)

! compute step size for angles and number of Legendre moments to compute based on Wiscombe's criteria
step = 2 / (numang - 1)
nmom = int(2 * (xx_dp + 4 * xx_dp**one_third + 2))
xx = real(xx_dp)
momdim = nmom + 1

! allocate arrays for MIEV0 calculation
allocate(xmu(numang))
allocate(s1(numang))
allocate(s2(numang))
allocate(pmom(0:momdim, 1))

! compute cosine of angles for scattering amplitude functions
do n = 1, numang
   xmu(n) = real(1 - (n - 1)*step)
end do

! perform Mie calculation
call MIEV0(xx, crefin, perfct, mimcut, anyang,       &
           numang, xmu, nmom, ipolzn, momdim, prnt,  &
           qext, qsca, gqsc, pmom, sforw, sback, s1, &
           s2, tforw, tback, spike)

! compute normalization factor for Legendre moments
fnorm  = 4._real64 / (xx_dp**2 * real(qsca, kind=real64))

! serialize results as json
write(*, fmt="(a)") "{"
write(*, fmt="(a,es15.8,a)") '  "refractive_index_real": ', mr, ','
write(*, fmt="(a,es15.8,a)") '  "refractive_index_imaginary": ', mi, ','
write(*, fmt="(a,es15.8,a)") '  "particle_radius": ', rad, ','
write(*, fmt="(a,es15.8,a)") '  "wavelength": ', lam, ','
write(*, fmt="(a,es15.8,a)") '  "size_parameter": ', xx, ','
write(*, fmt="(a,es15.8,a)") '  "extinction_efficiency": ', qext, ','
write(*, fmt="(a,es15.8,a)") '  "scattering_efficiency": ', qsca, ','
write(*, fmt="(a,es15.8,a)") '  "absorption_efficiency": ', qext - qsca, ','
write(*, fmt="(a,es15.8,a)") '  "sphere_geometric_cross_section": ', pi * real(rad, kind=real64)**2, ','
write(*, fmt="(a,es15.8,a)") '  "single_scattering_albedo": ', qsca/qext, ','
write(*, fmt="(a,es15.8,a)") '  "asymmetry_factor": ', gqsc/qsca, ','
write(*, fmt="(a,es15.8,a)") '  "radiation_pressure_efficiency": ', qext - gqsc, ','
write(*, fmt="(a,i4,a)") '  "number_of_moments": ', nmom, ','
write(*, fmt="(a)") '  "legendre_moments": ['
do i = 0, nmom-1
    write(*, fmt="(a,es15.8,a)") '    ', fnorm * real(pmom(i,1), kind=real64), ','
end do
write(*, fmt="(a,es15.8,a)") '    ', fnorm * real(pmom(nmom,1), kind=real64)
write(*, fmt="(a)") '  ]'
write(*, fmt="(a)") '}'

! free memory allocated arrays
deallocate(xmu)
deallocate(s1)
deallocate(s2)
deallocate(pmom)

end program miescat
