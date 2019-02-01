!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Written by:
!!
!!   Luis 'Tony' Martinez <tony.mtos@gmail.com> (Johns Hopkins University)
!!
!!   Copyright (C) 2012-2013, Johns Hopkins University
!!
!!   This file is part of The Actuator Turbine Model Library.
!!
!!   The Actuator Turbine Model is free software: you can redistribute it
!!   and/or modify it under the terms of the GNU General Public License as
!!   published by the Free Software Foundation, either version 3 of the
!!   License, or (at your option) any later version.
!!
!!   The Actuator Turbine Model is distributed in the hope that it will be
!!   useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU General Public License for more details.
!!
!!   You should have received a copy of the GNU General Public License
!!   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*******************************************************************************
module atm_input_util
!*******************************************************************************
! This module reads the input files for the actuator turbine model module (ATM)

! Module for dynamic allocation variables
use atm_base

implicit none

! The variables for the ATM are defined here
integer :: numberOfTurbines
integer :: outputInterval
integer :: updateInterval

! This type will store the necessary variables for -each- turbine
! To declare: type(turbineArray_t), allocatable, dimension(:) :: turbineArray
! To access the variables do: turbineArray(n) % variableNeeded
type turbineArray_t
    ! Variables that are read from input files
    character(128) :: turbineName ! Name of turbine ('turbine1')
    character(128) :: turbineType ! Name of turbine type ('NREL5MWRef')
    real(rprec), dimension(3) :: baseLocation ! Location of the base (0 0 0)
    integer :: numBladePoints ! Number of points along each blade
    character(128) :: bladeUpdateType !
    real(rprec) :: epsilon ! Width of the smearing Gaussian function
    character(128) :: sampling ! Sampling method for velocity atPoint or Spalart
    character(128) :: rotationDir ! Direction of rotation ('cw')
    real(rprec) :: Azimuth   ! Angle of rotation of the rotor
    real(rprec) :: RotSpeed  ! Speed of the rotor (rpm)
    real(rprec) :: Pitch
    real(rprec) :: NacYaw    ! The yaw angle of the nacelle
    real(rprec) :: fluidDensity   ! The density of the fluid (used for power)
    integer :: numAnnulusSections ! Number of annulus sections on each blade
    real(rprec) :: AnnulusSectionAngle ! Number of annulus sections on each blade
    real(rprec) :: deltaNacYaw = 0._rprec ! Change in nacelle angle
    real(rprec) :: TSR = 0._rprec ! Tip speed ratio
    real(rprec) :: PitchControlAngle = 0._rprec
    real(rprec) :: IntSpeedError = 0._rprec
    real(rprec) :: IntPowerError = 0._rprec
    logical :: tipALMCorrection = .false. ! Includes a correction for tip
    logical :: rootALMCorrection = .false. ! Includes a correction for tip
    real(rprec) :: optimalEpsilonChord = 0.25 ! The optimal epsilon / chord

    ! Not read variables
    real(rprec) :: thrust ! Total turbine thrust
    real(rprec) :: torqueRotor ! Rotor torque
    real(rprec) :: torqueGen ! Generator torque
    real(rprec) :: powerRotor ! Rotor Power
    real(rprec) :: powerGen ! Generator Power
    logical :: nacelle  ! Includes a nacelle yes or no
    real(rprec) :: nacelleEpsilon ! Width of the smearing Gaussian function
    real(rprec) :: nacelleCd = 0._rprec ! Drag coefficient for the nacelle
    real(rprec) :: VelNacelle_sampled = 0._rprec ! Sampled nacelle Vel
    real(rprec) :: VelNacelle_corrected = 0._rprec ! Corrected nacelle Vel
    real(rprec) :: u_infinity_mean = 0._rprec ! Mean velocity

    ! The MPI communicator for this turbine
    integer :: TURBINE_COMM_WORLD

    ! Master processor for turbine (will write output)
    integer :: master

    ! Flag to know if this turbine operates or not in this processor
    logical :: operate

    integer :: turbineTypeID ! Identifies the type of turbine

    !!-- Important geometry data
    ! Collection of all the actuator points (blade, annular section, point, 3)
    real(rprec), allocatable, dimension(:,:,:,:) :: bladePoints
    ! The solidity at each actuator section
    real(rprec), allocatable, dimension(:,:,:) :: solidity
    ! Collection of radius of each point (different because of coning)
    real(rprec), allocatable, dimension(:,:,:) :: bladeRadius
    ! Twist angle along the blade
    real(rprec), allocatable, dimension(:,:,:) :: twistAng
    ! Chord along the blade
    real(rprec), allocatable, dimension(:,:,:) :: chord
    ! Section type along the blade
    integer,     allocatable, dimension(:,:,:) :: sectionType
    ! Forces on each actuator point (blade, annular section, point, 3)
    real(rprec), allocatable, dimension(:,:,:,:) :: bladeForces
    ! Drag force of Nacelle
    real(rprec), dimension(3) :: nacelleForce
    ! Forces on each actuator point (blade, annular section, point, 3)
    real(rprec), allocatable, dimension(:,:,:,:) :: integratedBladeForces
    ! Vectors at each actuator point defining the local reference frame
    ! (blade, annular section, point, 3, 3) (three vectors)
    real(rprec), allocatable, dimension(:,:,:,:,:) :: bladeAlignedVectors
    ! The wind U projected onto the bladeAlignedVectors plus rotational speed
    ! (blade, annular section, point, 3, 3) (three vectors)
    real(rprec), allocatable, dimension(:,:,:,:) :: windVectors
    ! Angle of attack at each each actuator point
    real(rprec), allocatable, dimension(:,:,:) :: alpha
    ! Velocity magnitud at each each actuator point
    real(rprec), allocatable, dimension(:,:,:) :: Vmag
    ! Lift coefficient at each actuator point
    real(rprec), allocatable, dimension(:,:,:) :: Cl
    ! Lift coefficient correction at each actuator point
    real(rprec), allocatable, dimension(:,:,:) :: dalpha
    ! Drag coeficient at each actuator point
    real(rprec), allocatable, dimension(:,:,:) :: Cd
    ! Lift at each actuator point
    real(rprec), allocatable, dimension(:,:,:) :: lift
    ! Drag at each actuator point
    real(rprec), allocatable, dimension(:,:,:) :: drag
    ! Axial force at each actuator point
    real(rprec), allocatable, dimension(:,:,:) :: axialForce
    ! Tangential force at each actuator point
    real(rprec), allocatable, dimension(:,:,:) :: tangentialForce
    ! The function G=1/2 * Cl * c * Vmag^2
    real(rprec), allocatable, dimension(:,:,:) :: G
    ! The derivative of G
    real(rprec), allocatable, dimension(:,:,:) :: dG
    ! Cl base for the correction
    real(rprec), allocatable, dimension(:,:,:) :: Cl_b
    ! Slope of the lift curve base for the correction (dCl/ d alpha)
    real(rprec), allocatable, dimension(:,:,:) :: dCldalpha
    ! The induced velocity in the LES
    real(rprec), allocatable, dimension(:,:,:) :: uy_LES
    ! The induced velocity from the optimal solution
    real(rprec), allocatable, dimension(:,:,:) :: uy_opt
    ! The induced velocity vector from the LES
    real(rprec), allocatable, dimension(:,:,:,:) :: uy_LES_vec
    ! The induced drag velocity vector from the optimal
    real(rprec), allocatable, dimension(:,:,:,:) :: ux_LES_vec
    ! The induced velocity vector from the optimal
    real(rprec), allocatable, dimension(:,:,:,:) :: uy_opt_vec
    ! The inflow velocity in the relative frame of reference
    real(rprec), allocatable, dimension(:,:,:,:) :: Uinf_vec
    ! The change in induced velocity between epsilon_LES and epsilon_opt
    real(rprec), allocatable, dimension(:,:,:,:) :: du

    ! These variables are to make corrections based on optimum value of epsilon
    real(rprec), allocatable, dimension(:,:,:) :: epsilon_opt

!~     ! The circulation needed for the correction
!~     real(rprec), allocatable, dimension(:,:,:) :: Gamma

    ! Induction factor and u infinity
    real(rprec), allocatable, dimension(:,:,:) :: induction_a
    real(rprec), allocatable, dimension(:,:,:) :: u_infinity

    ! These are dummies meant to be used for parallelization
    ! bladeVectorDummy store quantities along the blades which are vectors
    ! such as Force
    ! bladeScalarDummy will store scalar quantities along the blades such
    ! as lift coefficitent, angle of attack, etc
    real(rprec), allocatable, dimension(:,:,:,:) :: bladeVectorDummy
    real(rprec), allocatable, dimension(:,:,:) :: bladeScalarDummy

    ! An indicator of shaft direction.  The convention is that when viewed
    ! from upwind, the rotor turns clockwise for positive rotation angles,
    ! regardless of if it is an upwind or downwind turbine.  uvShaft is
    ! found by subtracting the rotor apex location from the tower shaft
    ! intersection point.  This vector switches direciton depending on
    ! if the turbine is upwind or downwind, so this uvShaftDir multiplier
    ! makes the vector consistent no matter what kind of turbine
     real(rprec) :: uvShaftDir

    ! Define the vector along the shaft pointing in the direction of the wind
    real(rprec), dimension(3) :: uvShaft

    ! List of locations of the rotor apex relative to the origin (m)
    real(rprec), dimension(3) :: rotorApex

    ! List of locations of the intersection of the tower axis and the shaft
    ! centerline relative to the origin (m).
    real(rprec), dimension(3) :: towerShaftIntersect

     ! Unit vector pointing along the tower (axis of yaw).
    real(rprec), dimension(3) :: uvTower

    ! Width of the actuator section
    real(rprec), allocatable, dimension(:) :: db

    ! Sphere radius which defines a sphere from the center of the rotor and
    ! identifies the volume onto which forces are applied
    real(rprec) :: projectionRadius
    real(rprec) :: projectionRadiusNacelle
    real(rprec) :: sphereRadius

    ! Location of Nacelle
    real(rprec), dimension(3) :: nacelleLocation

    ! Change of azimuth angle every time-step
    real(rprec) :: deltaAzimuth=0.

end type turbineArray_t

! This type stores all the airfoils and their AOA, Cd, Cl values
type airfoilType_t
    character(128) :: airfoilName          ! The type of Airfoil ('Cylinder1')
    integer :: n                           ! Number of data points
    ! The maximum number of points is chosen to be 150. If airfoil data has
    ! more than this then this number should be modified!
    real(rprec), dimension(150) :: AOA    ! Angle of Attack
    real(rprec), dimension(150) :: Cd     ! Drag coefficient
    real(rprec), dimension(150) :: Cl     ! Lift coefficient
    real(rprec), dimension(150) :: Cm     ! Moment coefficient

end type airfoilType_t

type turbineModel_t
    character(128) :: turbineType ! The type of turbine ('NREL5MWRef')
    integer :: NumBl          ! Number of turbine blades
    integer :: NumSec         ! Number of sections
    real(rprec) :: TipRad ! Radius from the root to the tip of the blade
    real(rprec) :: HubRad ! Radius of the hub
    real(rprec) :: UndSling ! Undersling length [distance from teeter pin
                            ! to the rotor apex]
    real(rprec) :: OverHang ! Distance from yaw axis to rotor apex
    real(rprec) :: TowerHt  ! Tower height
    real(rprec) :: Twr2Shft ! Vertical distance from tower-top to rotor shaft
    real(rprec) :: ShftTilt ! Rotor shaft tilt angle
    real(rprec) :: PreCone  ! Angle which the blades are coned
    real(rprec) :: GBRatio  ! Gearbox Ratio
    real(rprec) :: GenIner  ! Generator inertia
    real(rprec) :: HubIner  ! Inertia of the hub
    real(rprec) :: BladeIner ! Inertia of the blades
    real(rprec) :: DriveTrainIner ! Inertia of the Drive train

    ! Blade section quantities (maximum number of sections 100, easy modify)
    real(rprec), dimension(100) :: chord, twist, radius
    integer, dimension(100) :: sectionType

    ! The airfoil type properties ( includes AOA, Cl, and Cd) Attempt 1
    type(airfoilType_t), allocatable, dimension(:) :: airfoilType

    ! Torque controller variables
    character(64) :: TorqueControllerType
    real(rprec) :: CutInGenSpeed
    real(rprec) :: RatedGenSpeed
    real(rprec) :: Region2StartGenSpeed
    real(rprec) :: Region2EndGenSpeed
    real(rprec) :: CutInGenTorque
    real(rprec) :: RatedGenTorque
    real(rprec) :: RateLimitGenTorque
    real(rprec) :: KGen
    real(rprec) :: TorqueControllerRelax

    ! Pitch controller variables
    character(64) :: PitchControllerType
    real(rprec) :: PitchControlAngleMax   ! Maximum pitch angle
    real(rprec) :: PitchControlAngleMin   ! Minimum pitch angle
    real(rprec) :: PitchControlAngleK     ! Angle at which sensitivity doubles
    real(rprec) :: PitchControlKP0        ! Proportional term at angle = 0
    real(rprec) :: PitchControlKI0        ! Integral term at angle = 0

    ! Yaw controller variables
    character(64) :: YawControllerType    ! Name of yaw controller type
    character(64) :: YawControllerFile    ! File that contains the yaw time info
    real(rprec), allocatable, dimension(:)  :: yaw_time  ! time for yaw (seconds)
    real(rprec), allocatable, dimension(:)  :: yaw_angle ! yaw angle (degrees)

end type turbineModel_t

! Declare turbine array variable
type(turbineArray_t), allocatable, dimension(:) , target :: turbineArray

! Declare turbine model variable (stores information for turbine models)
type(turbineModel_t), allocatable, dimension(:), target :: turbineModel

! Name of the utility used
character (*), parameter :: mod_name = 'atm_input_util'
character (*), parameter :: input_conf = './inputATM/turbineArrayProperties'
character (*), parameter :: comment = '!'
character (*), parameter :: block_entry = '{' ! The start of a block
character (*), parameter :: block_exit = '}' ! The end of a block
character (*), parameter :: equal = '='
character (*), parameter :: esyntax = 'syntax error at line'
character (*), parameter :: array_entry = '(' ! The start of an array
character (*), parameter :: array_exit  = ')' ! The end of an array

! Delimiters used for reading vectors and points
character(*), parameter :: delim_minor=','
character(*), parameter :: delim_major='//'

! Thresh hold for evaluating differences in floating point values.
!real(rprec), parameter :: thresh = 1.0e-6_rprec

! Variables used to read lines in file
integer :: block_entry_pos ! Determines if the line is the start of a block
integer :: block_exit_pos ! Determines if the line is the end of a block
integer :: array_entry_pos ! Determines if the line is the start of a block
integer :: array_exit_pos ! Determines if the line is the end of a block
integer :: equal_pos ! Determines if there is an equal sign
integer :: ios
logical :: exst ! Used to determine existence of a file

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_input_conf()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

character (*), parameter :: sub_name = mod_name // '.read_input_conf'
integer :: n = 0 ! Counter for the wind turbines
integer :: lun =1 ! Reference number for input file
integer :: line ! Counts the current line in a file
character (128) :: buff ! Stored the read line

! Check that the configuration file exists
inquire (file=input_conf, exist=exst)

! Open file
if (exst) then
    ! Open the input file
    open (lun, file=input_conf, action='read')
else
    ! Error for non existing file
    call error ('file ' // input_conf // ' does not exist')
endif

! Read the file line by line *Counter starts at 0 and modified inside subroutine
line = 0
do
! Read line by line (lun=file number)
    call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                   array_entry_pos, array_exit_pos, equal_pos, ios )

    if (ios /= 0) exit ! Exit if reached end of file

    ! This will read the numberOfTurbines integer
        if( buff(1:16) == 'numberOfTurbines' ) then
            read(buff(17:), *) numberOfTurbines
!            write(*,*) 'numberOfTurbines is: ', numberOfTurbines
            ! Allocate space for the wind turbine variables
            allocate(turbineArray(numberOfTurbines))
            cycle
        else if( buff(1:14) == 'outputInterval' ) then
            read(buff(15:), *) outputInterval
!            write(*,*)  'outputInterval is: ', outputInterval
            cycle
        else if( buff(1:14) == 'updateInterval' ) then
            read(buff(15:), *) updateInterval
!            write(*,*)  'updateInterval is: ', updateInterval
            cycle
        endif

    if (block_entry_pos /= 0) then ! This will start reading turbine block
        n = n + 1     ! Increment turbine counter
        if (n .gt. numberOfTurbines) exit
        ! Read the name of the turbine
        read(buff(1:index(buff, block_entry)-1), *) turbineArray(n) &
        % turbineName
    endif
    if (block_entry_pos == 0) then ! This will start reading turbine block
        if( buff(1:11) == 'turbineType' ) then
            read(buff(12:), *) turbineArray(n) % turbineType
!            write(*,*)  'turbineType is: ', turbineArray(n) % turbineType
        endif
        if( buff(1:12) == 'baseLocation' ) then
            read(buff(13:), *) turbineArray(n) % baseLocation
!            write(*,*)  'baseLocation is: ', turbineArray(n) % baseLocation
        endif
        if( buff(1:14) == 'numBladePoints' ) then
            read(buff(15:), *) turbineArray(n) % numBladePoints
!            write(*,*)  'numBladePoints is: ', turbineArray(n) % numBladePoints
            ! Allocation depending on the number of blade points
        endif
        if( buff(1:7) == 'epsilon' ) then
            read(buff(8:), *) turbineArray(n) % epsilon
!            write(*,*)  'epsilon is: ', turbineArray(n) % epsilon
        endif
        if( buff(1:8) == 'sampling' ) then
            read(buff(9:), *) turbineArray(n) % sampling
!~             write(*,*)  'sampling is: ', turbineArray(n) % sampling
        endif
        if( buff(1:11) == 'rotationDir' ) then
            read(buff(12:), *) turbineArray(n) % rotationDir
!            write(*,*)  'rotationDir is: ', turbineArray(n) % rotationDir
        endif
        if( buff(1:7) == 'Azimuth' ) then
            read(buff(8:), *) turbineArray(n) % Azimuth
!            write(*,*)  'Azimuth is: ', turbineArray(n) % Azimuth
        endif
        if( buff(1:8) == 'RotSpeed' ) then
            read(buff(9:), *) turbineArray(n) % RotSpeed
!            write(*,*)  'RotSpeed is: ', turbineArray(n) % RotSpeed
        endif
        if( buff(1:5) == 'Pitch' ) then
            read(buff(6:), *) turbineArray(n) % Pitch
!            write(*,*)  'Pitch is: ', turbineArray(n) % Pitch
        endif
        if( buff(1:6) == 'NacYaw' ) then
            read(buff(7:), *) turbineArray(n) % NacYaw
!            write(*,*)  'NacYaw is: ', turbineArray(n) % NacYaw
        endif
        if( buff(1:12) == 'fluidDensity' ) then
            read(buff(13:), *) turbineArray(n) % fluidDensity
!            write(*,*)  'fluidDensity is: ', turbineArray(n) % fluidDensity
        endif
        if( buff(1:18) == 'numAnnulusSections' ) then
            read(buff(19:), *) turbineArray(n) % numAnnulusSections
!            write(*,*)  'numAnnulusSections is: ', &
!                         turbineArray(n) % numAnnulusSections
        endif
        if( buff(1:19) == 'annulusSectionAngle' ) then
            read(buff(20:), *) turbineArray(n) % annulusSectionAngle
!            write(*,*)  'annulusSectionAngle is: ', &
!                         turbineArray(n) % annulusSectionAngle
        endif
        if( buff(1:11) == 'nacelleFlag' ) then
            read(buff(12:), *) turbineArray(n) % nacelle
!            write(*,*)  'annulusSectionAngle is: ', &
!                         turbineArray(n) % annulusSectionAngle
        endif
        if( buff(1:9) == 'nacelleCd' ) then
            read(buff(10:), *) turbineArray(n) % nacelleCd
!~             write(*,*)  'cd is: ', &
!~                          turbineArray(n) % nacelleCd
        endif
        if( buff(1:14) == 'nacelleEpsilon' ) then
            read(buff(15:), *) turbineArray(n) % nacelleEpsilon
!~             write(*,*)  'cd is: ', &
!~                          turbineArray(n) % nacelleEpsilon
        endif
        if( buff(1:3) == 'TSR' ) then
            read(buff(4:), *) turbineArray(n) % TSR
!~             write(*,*)  'TSR is: ', &
!~                          turbineArray(n) % TSR
        endif
        if( buff(1:16) == 'tipALMCorrection' ) then
            read(buff(17:), *) turbineArray(n) % tipALMCorrection
!~             write(*,*)  'tipALMCorrection is: ', &
!~                          turbineArray(n) % tipALMCorrection
        endif
        if( buff(1:17) == 'rootALMCorrection' ) then
            read(buff(18:), *) turbineArray(n) % rootALMCorrection
!~             write(*,*)  'rootALMCorrection is: ', &
!~                          turbineArray(n) % rootALMCorrection
        endif
        if( buff(1:19) == 'optimalEpsilonChord' ) then
            read(buff(20:), *) turbineArray(n) % optimalEpsilonChord
!~             write(*,*)  'optimalEpsilonChord is: ', &
!~                          turbineArray(n) % optimalEpsilonChord
        endif
    endif
end do

if( .not. allocated( turbineArray ) ) then
    write(*,*) 'Did not allocate memory for turbineArray'
stop
endif
close (lun)

call read_turbine_model_variables ()

end subroutine read_input_conf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_turbine_model_variables ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
implicit none
integer :: i, j, c ! counter
integer :: numTurbinesDistinct ! Number of different turbine types
character(128) :: currentTurbineType ! Will store turbineType in loop
character(139) :: input_turbine
integer :: lun =19  ! Reference number for input file
integer :: line ! Counts the current line in a file
character (128) :: buff ! Stored the read line
integer:: numAirfoils ! Number of distinct airfoils
integer :: k, p ! Used to loop through aifoil types and character counter
integer :: numAnnulusSections, numBladePoints, numBl, numSec

! Variables for reading yaw control file
integer :: ios, N
integer :: yawfile = 29

! Name of all the airfoils types (max 20) If more than this increase the number
character(128), dimension(20) :: airfoils

! Initialize variables for the loop
! Will find the number of distincet turbines to allocate memory for turbineModel
numTurbinesDistinct=0
currentTurbineType=turbineArray(1) % turbineType

! Counter for number of distinct turbines
c = 0
! Find how many turbine types there are
do i=1,numberOfTurbines

    ! Find if the name is repeated
    do j=1,i-1
        if (turbineArray(i) % turbineType .eq. turbineArray(j) % turbineType) then
            c = 1
            turbineArray(i) % turbineTypeID = turbineArray(j) % turbineTypeID
        endif
    enddo

    ! Assign a new turbine if c = 0
    if (c .eq. 0) then
        numTurbinesDistinct = numTurbinesDistinct + 1
        turbineArray(i) % turbineTypeID = numTurbinesDistinct
    endif

    ! Restart the counter at 0
    c = 0
!~     write(*,*) 'turbine ',i,' is type ',turbineArray(i) % turbineType, ' ', turbineArray(i) % turbineTypeID
enddo

! Allocate space for turbine model variables
allocate(turbineModel(numTurbinesDistinct))
write(*,*) 'Distinct Turbines=', numTurbinesDistinct

! This will store the turbine types on each turbine model ("NREL5MW")
!~ numTurbinesDistinct=1
!~ currentTurbineType=turbineArray(1) % turbineType
!~ turbineModel(numTurbinesDistinct) % turbineType = turbineArray(1) % turbineType
!~ do i=1,numberOfTurbines
!~ currentTurbineType=turbineArray(i) % turbineType ! added
!~     if (turbineArray(i) % turbineType .ne. currentTurbineType) then
!~     numTurbinesDistinct = numTurbinesDistinct+1
!~     turbineModel(numTurbinesDistinct) % turbineType = &
!~     turbineArray(i) % turbineType
!~     endif
!~ enddo
!~ numTurbinesDistinct=1
!~ currentTurbineType=turbineArray(1) % turbineType
!~ turbineModel(numTurbinesDistinct) % turbineType = turbineArray(1) % turbineType

do i=1,numberOfTurbines
    turbineModel(turbineArray(i) % turbineTypeID) % turbineType = &
    turbineArray(i) % turbineType
enddo

! Read the input properties for each turbine type
do i = 1, numTurbinesDistinct

    input_turbine = './inputATM/' // turbineModel(i) % turbineType

    ! Check that the configuration file exists
    inquire (file=input_turbine, exist=exst)

!~  write(*,*) 'Error NOT Here', i, lun, input_turbine, exst
    ! Open file
    if (exst) then
        ! Open the input file
        open (lun, file=input_turbine, action='read')
!~  write(*,*) 'Error Here'

    else
        ! Error for non existing file
        call error ('file ' // input_turbine // ' does not exist')
    endif

!    write(*,*) 'Reading Turbine Model Properties for: ', &
!                turbineModel(i) % turbineType
    ! Read the file line by line - starts at 0 and modified inside subroutine
    line = 0
    do
    ! Read line by line (lun=file number)
        call readline( lun, line, buff, block_entry_pos, block_exit_pos, &
                       array_entry_pos, array_exit_pos, equal_pos, ios )
        if (ios /= 0) exit
        ! This will all the input variables
        if( buff(1:5) == 'NumBl' ) then
            read(buff(6:), *) turbineModel(i) % NumBl
!            write(*,*) 'NumBl is: ', turbineModel(i) % NumBl
        endif
        if( buff(1:6) == 'TipRad' ) then
            read(buff(7:), *) turbineModel(i) % TipRad
!            write(*,*) 'TipRad is: ', turbineModel(i) % TipRad
        endif
        if( buff(1:6) == 'HubRad' ) then
            read(buff(7:), *) turbineModel(i) % HubRad
!            write(*,*) 'HubRad is: ', turbineModel(i) % HubRad
        endif
        if( buff(1:8) == 'UndSling' ) then
            read(buff(9:), *) turbineModel(i) % UndSling
!            write(*,*) 'UndSling is: ', turbineModel(i) % UndSling
        endif
        if( buff(1:8) == 'OverHang' ) then
            read(buff(9:), *) turbineModel(i) % OverHang
!            write(*,*) 'OverHang is: ', turbineModel(i) % OverHang
        endif
        if( buff(1:7) == 'TowerHt' ) then
            read(buff(8:), *) turbineModel(i) % TowerHt
!            write(*,*) 'TowerHt is: ', turbineModel(i) % TowerHt
        endif
        if( buff(1:8) == 'Twr2Shft' ) then
            read(buff(9:), *) turbineModel(i) % Twr2Shft
!            write(*,*) 'Twr2Shft is: ', turbineModel(i) % Twr2Shft
        endif
        if( buff(1:8) == 'ShftTilt' ) then
            read(buff(9:), *) turbineModel(i) % ShftTilt
!            write(*,*) 'ShftTilt is: ', turbineModel(i) % ShftTilt
        endif
        if( buff(1:7) == 'PreCone' ) then
            read(buff(8:), *) turbineModel(i) % PreCone
!            write(*,*) 'PreCone is: ', turbineModel(i) % PreCone
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This will read the torque controller type and its properties
        if( buff(1:20) == 'TorqueControllerType' ) then
            read(buff(21:), *) turbineModel(i) % TorqueControllerType
!            write(*,*) 'TorqueControllerType is: ', turbineModel(i) % TorqueControllerType
        endif
        if( buff(1:13) == 'CutInGenSpeed' ) then
            read(buff(14:), *) turbineModel(i) % CutInGenSpeed
!            write(*,*) 'CutInGenSpeed is: ', turbineModel(i) % CutInGenSpeed
        endif
        if( buff(1:13) == 'RatedGenSpeed' ) then
            read(buff(14:), *) turbineModel(i) % RatedGenSpeed
!            write(*,*) 'RatedGenSpeed is: ', turbineModel(i) % RatedGenSpeed
        endif
        if( buff(1:20) == 'Region2StartGenSpeed' ) then
            read(buff(21:), *) turbineModel(i) % Region2StartGenSpeed
!            write(*,*) 'Region2StartGenSpeed is: ', turbineModel(i) % Region2StartGenSpeed
        endif
        if( buff(1:18) == 'Region2EndGenSpeed' ) then
            read(buff(19:), *) turbineModel(i) % Region2EndGenSpeed
!            write(*,*) 'Region2EndGenSpeed is: ', turbineModel(i) % Region2EndGenSpeed
        endif
        if( buff(1:14) == 'CutInGenTorque' ) then
            read(buff(15:), *) turbineModel(i) % CutInGenTorque
!            write(*,*) 'CutInGenTorque is: ', turbineModel(i) % CutInGenTorque
        endif
        if( buff(1:14) == 'RatedGenTorque' ) then
            read(buff(15:), *) turbineModel(i) % RatedGenTorque
!            write(*,*) 'RatedGenTorque is: ', turbineModel(i) % RatedGenTorque
        endif
        if( buff(1:18) == 'RateLimitGenTorque' ) then
            read(buff(19:), *) turbineModel(i) % RateLimitGenTorque
!            write(*,*) 'RateLimitGenTorque is: ', turbineModel(i) % RateLimitGenTorque
        endif
        if( buff(1:4) == 'KGen' ) then
            read(buff(5:), *) turbineModel(i) % KGen
!            write(*,*) 'KGen is: ', turbineModel(i) % KGen
        endif
        if( buff(1:21) == 'TorqueControllerRelax' ) then
            read(buff(22:), *) turbineModel(i) % TorqueControllerRelax
!            write(*,*) 'TorqueControllerRelax is: ', turbineModel(i) % TorqueControllerRelax
        endif
        if( buff(1:7) == 'GBRatio' ) then
            read(buff(8:), *) turbineModel(i) % GBRatio
!            write(*,*) 'GBRatio is: ', turbineModel(i) % GBRatio
        endif
        if( buff(1:9) == 'BladeIner' ) then
            read(buff(10:), *) turbineModel(i) % BladeIner
!            write(*,*) 'BladeIner is: ', turbineModel(i) % BladeIner
        endif
        if( buff(1:7) == 'HubIner' ) then
            read(buff(8:), *) turbineModel(i) % HubIner
!            write(*,*) 'HubIner is: ', turbineModel(i) % HubIner
        endif
        if( buff(1:7) == 'GenIner' ) then
            read(buff(8:), *) turbineModel(i) % GenIner
!            write(*,*) 'GenIner is: ', turbineModel(i) % GenIner
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This will read the pitch control values
        if( buff(1:19) == 'PitchControllerType' ) then
            read(buff(20:), *) turbineModel(i) % PitchControllerType
            write(*,*) 'PitchControllerType is: ', turbineModel(i) % PitchControllerType
        endif
        if( buff(1:20) == 'PitchControlAngleMax' ) then
            read(buff(21:), *) turbineModel(i) % PitchControlAngleMax
            write(*,*) 'PitchControlAngleMax is: ', turbineModel(i) % PitchControlAngleMax
        endif
        if( buff(1:20) == 'PitchControlAngleMin' ) then
            read(buff(21:), *) turbineModel(i) % PitchControlAngleMin
            write(*,*) 'PitchControlAngleMin is: ', turbineModel(i) % PitchControlAngleMin
        endif
        if( buff(1:18) == 'PitchControlAngleK' ) then
            read(buff(19:), *) turbineModel(i) % PitchControlAngleK
            write(*,*) 'PitchControlAngleK is: ', turbineModel(i) % PitchControlAngleK
        endif
        if( buff(1:15) == 'PitchControlKP0' ) then
            read(buff(16:), *) turbineModel(i) % PitchControlKP0
            write(*,*) 'PitchControlKP0 is: ', turbineModel(i) % PitchControlKP0
        endif
        if( buff(1:15) == 'PitchControlKI0' ) then
            read(buff(16:), *) turbineModel(i) % PitchControlKI0
            write(*,*) 'PitchControlKI0 is: ', turbineModel(i) % PitchControlKI0
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This will read the yaw control values
        if( buff(1:17) == 'YawControllerType' ) then
            read(buff(18:), *) turbineModel(i) % YawControllerType
            write(*,*) 'YawControllerType is: ', turbineModel(i) % YawControllerType
        endif
        if( buff(1:17) == 'YawControllerFile' ) then
            read(buff(18:), *) turbineModel(i) % YawControllerFile
            write(*,*) 'YawControllerFile is: ', turbineModel(i) % YawControllerFile
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This will read the airfoils
        if ( buff(1:8) == 'Airfoils' ) then ! Start reading airfoil block
            numAirfoils=0 ! Conuter for the number of distince airfoils
            array_entry_pos=0 ! If 'Airfoils(' then make array_entry_pos zero

            do while (array_entry_pos == 0)
                call readline( lun, line, buff, block_entry_pos,              &
                block_exit_pos, array_entry_pos, array_exit_pos, equal_pos, ios)
                if (ios /= 0) exit        ! exit if end of file reached

                if (array_entry_pos /= 0) then
                    call readline( lun, line, buff, block_entry_pos,          &
                                   block_exit_pos, array_entry_pos,           &
                                   array_exit_pos, equal_pos, ios)
                   if (ios /= 0) exit        ! exit if end of file reached
                endif

                if (array_exit_pos /= 0) exit     ! exit if end of file reached
                numAirfoils = numAirfoils + 1     ! Increment airfoil counter
                airfoils(numAirfoils)=buff ! Stores the name of the airfoil
            enddo

            ! Allocate the airfoilTypes
            allocate(turbineModel(i) % airfoilType(numAirfoils))

            ! Loop through all the airfoil types and look-up lift and drag
            do k=1,numAirfoils
                call eat_whitespace (airfoils(k)) ! Eliminate white space
                p=len(trim(airfoils(k)))-1        ! Length without last element
                ! Airfoil type (2:p) is used to eliminate the commas ""
                turbineModel(i) % airfoilType % airfoilName =  airfoils(k)(2:p)
                ! Read each airfoil accordingly
                call read_airfoil( turbineModel(i) % airfoilType(k) )
            enddo
        endif  ! End of Airfoil loop

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This will read the blade properties
        if ( buff(1:9) .eq. 'BladeData' ) then ! Start reading blade data block
                NumSec = 0
                do while (array_exit_pos .eq. 0)
                    read (lun, '(a)', iostat=ios) buff ! Read the comment line
                    if (scan(buff,'(' ) .ne. 0) cycle
                    if (scan(buff,'!' ) .ne. 0) cycle
                    if (scan(buff,')' ) .ne. 0) exit

                    ! Number of sections will account for all blade sections
                    NumSec=NumSec+1
                    turbineModel(i) % NumSec = NumSec

                    ! Read in radius, chord, twist, type
                    read(buff,*) turbineModel(i) % radius(NumSec),             &
                    turbineModel(i) % chord(NumSec),   &
                    turbineModel(i) % twist(NumSec),   &
                    turbineModel(i) % sectionType(NumSec)

                    ! Add one to airfoil identifier. List starts at 0, now
                    ! it will start at 1
                    turbineModel(i) % sectionType( NumSec )  =                 &
                    turbineModel(i) % sectionType( NumSec ) + 1
                enddo
        endif
    enddo
    close (lun)

    ! Read the yaw control file, if applicable
    if ( turbineModel(i) % YawControllerType == "timeYawTable" ) then
        open (unit = yawfile, file="inputATM/"//                               &
                trim(turbineModel(i) % YawControllerFile),                     &
                form = "formatted", status = "old", action = "read")

        ! Determine number of lines
        ios = 0
        N = -1
        do while (ios == 0)
            N = N + 1
            read(unit = yawfile, fmt = *, iostat = ios)
        enddo
        rewind(unit = yawfile)

        ! Allocate variables
        allocate(turbineModel(i) % yaw_time(N))
        allocate(turbineModel(i) % yaw_angle(N))

        ! now read the variables and close file
        do j = 1, N
            read(unit = yawfile, fmt = *) turbineModel(i) % yaw_time(j),       &
                    turbineModel(i) % yaw_angle(j)
            print *, turbineModel(i) % yaw_time(j), &
                     turbineModel(i) % yaw_angle(j)
        enddo
        close(unit = yawfile)
    endif

    ! Calculate drive train inertia
    turbineModel(i) % DriveTrainIner = (real(turbineModel(i) % NumBl,rprec)) * &
              (turbineModel(i) % BladeIner) + (turbineModel(i) % HubIner) +    &
              ( turbineModel(i) % GBRatio ) * ( turbineModel(i) % GBRatio) *   &
              ( turbineModel(i) % GenIner )

enddo

! Allocate variables inside turbineArray
do  i=1,numberOfTurbines
numBladePoints = turbineArray(i) % numBladePoints
numAnnulusSections = turbineArray(i) % numAnnulusSections
j=turbineArray(i) % turbineTypeID
numBl=turbineModel(j) % numBl

    allocate(turbineArray(i) % bladeForces(numBl,                              &
             numAnnulusSections, numBladePoints,3) )
    allocate(turbineArray(i) % integratedBladeForces(numBl,                    &
             numAnnulusSections, numBladePoints,3) )
    allocate(turbineArray(i) % bladeAlignedVectors(numBl,                      &
             numAnnulusSections, numBladePoints,3,3) )
    allocate(turbineArray(i) % windVectors(numBl,                              &
             numAnnulusSections, numBladePoints,3) )
    allocate(turbineArray(i) % alpha(numBl,                                    &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % Vmag(numBl,                                     &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % Cl(numBl,                                       &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % dalpha(numBl,                                   &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % Cd(numBl,                                       &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % lift(numBl,                                     &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % drag(numBl,                                     &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % axialForce(numBl,                               &
             numAnnulusSections, numBladePoints))
    allocate(turbineArray(i) % tangentialForce(numBl,                          &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % induction_a(numBl,                              &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % u_infinity(numBl,                               &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % chord(numBl,                                    &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % twistAng(numBl,                                 &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % sectionType(numBl,                              &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % epsilon_opt(numBl,                              &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % uy_LES(numBl,                                   &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % uy_opt(numBl,                                   &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % Uinf_vec(numBl,                                 &
             numAnnulusSections, numBladePoints,3) )
    allocate(turbineArray(i) % uy_LES_vec(numBl,                               &
             numAnnulusSections, numBladePoints,3) )
    allocate(turbineArray(i) % ux_LES_vec(numBl,                               &
             numAnnulusSections, numBladePoints,3) )
    allocate(turbineArray(i) % uy_opt_vec(numBl,                               &
             numAnnulusSections, numBladePoints,3) )
    allocate(turbineArray(i) % du(numBl,                                       &
             numAnnulusSections, numBladePoints,3) )
    allocate(turbineArray(i) % G(numBl,                                        &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % dG(numBl,                                       &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % Cl_b(numBl,                                     &
             numAnnulusSections, numBladePoints) )
    allocate(turbineArray(i) % dCldalpha(numBl,                                &
             numAnnulusSections, numBladePoints) )

    ! Variables meant for parallelization
    allocate(turbineArray(i) % bladeVectorDummy(numBl,                         &
             numAnnulusSections, numBladePoints,3) )
    allocate(turbineArray(i) % bladeScalarDummy(numBl,                         &
             numAnnulusSections, numBladePoints) )

    ! Initialize to zero
    turbineArray(i) % du(:,:,:,:) = 0._rprec
    turbineArray(i) % Uinf_vec(:,:,:,:) = 0._rprec

enddo

end subroutine read_turbine_model_variables

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_print_initialize( )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Outputs the initialization to the screen
write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
if (numberOfTurbines==1) then
write(*,*) 'Actuator Turbine Model has been implemented with 1 turbine'
else
write(*,*) 'Actuator Turbine Model has been implemented with',  &
            numberOfTurbines ,'turbines'
endif
write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

end subroutine atm_print_initialize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_airfoil( airfoilType )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine reads the angle of attack, lift and drag for a specific
! airfoil type
integer :: lun=17 ! File identifier
integer :: q ! Counter
character(128) :: input_airfoil
type(airfoilType_t), intent(inout) :: airfoilType
real(rprec) :: AOA, Cd, Cl, Cm
! Name of the input file to read
input_airfoil= './inputATM/AeroData/' //trim (airfoilType % airfoilName)  &
               // '.dat'
!write (*,*) input_airfoil
! Open airfoil input file
open (lun, file=input_airfoil, action='read')

! Skip the first 14 lines of the Fortran input file
do q=1,14
    read(lun,*)
enddo

q=0 ! Initialize the counter back to 1 for the first element in the list

AOA=-181.
do while (AOA .lt. 180.00)
    q=q+1
    read(lun,*) AOA, Cl , Cd, Cm
    airfoilType % AOA(q) = AOA
    airfoilType % Cd(q) = Cd
    airfoilType % Cl(q) = Cl
    airfoilType % Cm(q) = Cm
enddo
airfoilType % n = q

close(lun)

end subroutine read_airfoil

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine readline(lun, line, buff, block_entry_pos, block_exit_pos, &
                    array_entry_pos, array_exit_pos, equal_pos, ios )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This subroutine reads the specified line and determines the attributes
! of the contents of the line.
!
implicit none

integer, intent(in) :: lun
integer, intent(inout) :: line

character(*), intent(inout) :: buff
integer, intent(out) :: block_entry_pos, block_exit_pos, equal_pos, ios, &
                        array_entry_pos, array_exit_pos

block_entry_pos = 0
block_exit_pos = 0
equal_pos = 0
ios = -1

do
    line = line + 1
    read (lun, '(a)', iostat=ios) buff

    if (ios /= 0) exit

    call eat_whitespace (buff)

    if (verify (buff, ' ') == 0) cycle  !--drop blank lines

    if (buff (1:len (comment)) == comment) cycle  !--drop comment lines

    block_entry_pos = index( buff, block_entry )
    block_exit_pos  = index( buff, block_exit )
    array_entry_pos = index( buff, array_entry )
    array_exit_pos  = index( buff, array_exit )
    equal_pos       = index( buff, equal )
    exit
enddo

return
end subroutine readline

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine eat_whitespace (buff, whtspc)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! eats leading and intermediate whitespace, fill trailing space with
! blanks
!
implicit none

character (*), intent (inout) :: buff
character (*), intent (in), optional :: whtspc
character (*), parameter :: whtspc_default = achar (9) // achar (32)
                            !--add more characters here if needed
character (1), parameter :: fill_char = ' '
character (1) :: tmp (len (buff))
character (1) :: fill (len (buff))

fill = fill_char
tmp = transfer (buff, tmp)

if (present (whtspc)) then
  tmp = pack (tmp, scan (tmp, whtspc) == 0, fill)
else
  tmp = pack (tmp, scan (tmp, whtspc_default) == 0, fill)
end if

buff = transfer (tmp, buff)

end subroutine eat_whitespace

end module atm_input_util
