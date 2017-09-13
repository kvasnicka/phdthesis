module mod_tools
!This module contains tools which may be useful in various projects.

use mod_types
use mod_IO
use IFPORT !Intel portability routine

implicit none

contains

!Function time_real_acc returns unix time with added miliseconds
!as a real number. It is used for tracking runtime as cpu_time and related
!functions don't work properly when multiple cores are used.
real(dp) function time_real_acc()
        integer(long_int_kind) :: time_int
        
        character(8) :: date_dat
        character(10) :: time_dat
        character(10) :: zone_dat
        integer, dimension(8) :: values_dat
        
        !call system_clock(count=time_int)
        !Even though this is an intrinsic subroutine so it would be nice to use it,
        !for some reason system_clock returned very inaccurate results (fluctuating by hundreds of seconds)
        !The portability function from IFORT library is used.
        time_int = time()

        call date_and_time(date_dat,time_dat,zone_dat,values_dat)
              
        time_real_acc = time_int + values_dat(8)/1000.0_dp
        
end function time_real_acc

!________________________________________________________________________________
!subroutine runtime_report

!This subroutine returns elapsed time as character(256).
!
!If convert_time == .true., this will be converted to days, hours, ..., otherwise
!it will be kept in seconds.


subroutine runtime_report(initial_time,convert_time,elapsed_time_string)

        implicit none

        real(dp), intent(in) :: initial_time
        logical, intent(in) :: convert_time
        character(256), intent(out) :: elapsed_time_string
        
        real(dp) :: cpu_time_elapsed = 0.0_dp, seconds = 0.0_dp, remainder = 0.0_dp
        integer :: days=0, hours=0, minutes=0
        character(256) :: tmp_string
                    
        
        cpu_time_elapsed = time_real_acc() - initial_time
         
        !write(*,*) cpu_time_elapsed !debug

        !Generate a string for output (gradually ammend things)
        elapsed_time_string = 'Elapsed time :'
        
        if (convert_time == .false.) then
                seconds = cpu_time_elapsed
                go to 101
        end if

        !This is not an efficient algorithm but it doesn't matter as 
        !this subroutine is not used very often
        days = floor(cpu_time_elapsed) / 86400 !(integer division)
        remainder = cpu_time_elapsed - days*86400
        
        hours = floor(remainder) / 3600
        remainder = remainder - hours*3600                
        
        minutes = floor(remainder) / 60
        seconds = remainder - minutes*60
        
        
        
        
        if(days>0) then
                write(tmp_string,'(I5)') days !assumed less than 99,999 days (more than my remaining lifetime so it should be sufficient)
                elapsed_time_string = trim(elapsed_time_string)//' '//adjustl(tmp_string)
                
                if(days>1) then
                        tmp_string = 'days'
                else
                        tmp_string = 'day'
                end if
                elapsed_time_string = trim(elapsed_time_string)//' '//adjustl(tmp_string)
        end if
        
        if(hours>0) then
                write(tmp_string,'(I2)') hours
                elapsed_time_string = trim(elapsed_time_string)//' '//adjustl(tmp_string)
                
                if(hours>1) then
                        tmp_string = 'hours'
                else
                        tmp_string = 'hour'
                end if
                elapsed_time_string = trim(elapsed_time_string)//' '//adjustl(tmp_string)
        end if
        
        if(minutes>0) then
                write(tmp_string,'(I2)') minutes
                elapsed_time_string = trim(elapsed_time_string)//' '//adjustl(tmp_string)
                
                if(minutes>1) then
                        tmp_string = 'minutes'
                else
                        tmp_string = 'minute'
                end if
                elapsed_time_string = trim(elapsed_time_string)//' '//adjustl(tmp_string)
        end if        
        
 
        !If elapsed time larger than 1 minute then display 2 decimals, otherwise just the whole number (floor)
        101 if(minutes > 0 .or. hours > 0 .or. days > 0 .or. seconds/60.0_dp > 1) then
                !seconds/60.0_dp > 1 added for the case of convert_time == .false. where the other variables are kept at initial value (0)
              write(tmp_string,'(I10)') floor(seconds) !I10 because of possible high values when convert_time == .false. 
        else
                write(tmp_string,'(F12.3)') seconds
        end if
        
        elapsed_time_string = trim(elapsed_time_string)//' '//adjustl(tmp_string)
        
        if(seconds>1) then
                tmp_string = 'seconds'
        else
                tmp_string = 'second'
        end if
        elapsed_time_string = trim(elapsed_time_string)//' '//adjustl(tmp_string)

end subroutine runtime_report
!________________________________________________________________________________
FUNCTION markov(S_0,P,T)
        integer :: S_0, T, seed
        real(kind=dp), dimension(:,:) :: P
        integer, dimension(T) :: markov

        integer, dimension(2):: dim_of_P !dimension of transition matrix (rows, columns)

        integer:: i,j,k !indices used in loops

        real(kind=dp):: sum_of_row

        real(kind=dp), dimension(:,:), ALLOCATABLE :: D

        real(kind=dp), dimension(T) :: shocks_under !underlying shocks from U[0,1]

        !First we want to check the inputs.
        !Check that the transition matrix P is a square matrix
        dim_of_P(1) = size(P,1)
        dim_of_P(2) = size(P,2)



        if(dim_of_P(1) /=dim_of_P(2)) then
               stop 'Error: The transition matrix P is not a square matrix.'
        endif



        !Check that all elements of P are non-negative
        do i = 1,dim_of_P(1)
                do j = 1,dim_of_P(2)
                        if (P(i,j) < 0) then
                                stop 'Error: Some element of transition matrix P is negative.'
                        end if
                end do
        end do


        !Check that all rows of the transition matrix sum to 1
        do i =1,dim_of_P(1) !cycling over rows
                sum_of_row = 0.0_dp !initialize sum at 0
                do j = 1,dim_of_P(2) !cycling over columns
                        sum_of_row = sum_of_row + P(i,j)
                end do
                if(sum_of_row /= 1.0_dp) then
                        stop 'Error: Rows of transition matrix must sum to 1.'
                end if
        end do



        !S_0 must be an index between 1 and the row/column dimension of transition matrix
        if (S_0 < 1 .or. S_0 > dim_of_P(1)) then
                stop 'Initial index S_0 outside of range.'
        end if

        !T must be positive
        if (T<0) then
                stop 'T (number of simulated realizations) must be positive.'
        end if
        !End of checking of inputs.


        !Division matrix D will have the same dimension as P (the first column will be the same)
        ALLOCATE(D(dim_of_P(1),dim_of_P(2)))
        do i=1,dim_of_P(1) !cycle over rows
                do j=1,dim_of_P(2) !cycle over columns
                        D(i,j) = 0_dp !initialize
                        do k = 1,j !add all elements up to column j to get c.d.f.
                                D(i,j) = D(i,j) + P(i,k)
                        end do
                end do
        end do



        !Now we have the division matrix, we can generate the shocks. First generate T realizations from U[0,1]
        do i=1,T
               CALL RANDOM_NUMBER(shocks_under)
        end do



        do i=1,T
                !first element treated separately (S_0 is not part of markov)
                if (i==1) then
                        do j=1,dim_of_P(2)
                                !index of initial shock is S_0 so we use that row of division matrix
                                if (shocks_under(i)<D(S_0,j)) then
                                        markov(i) = j
                                        exit
                                end if
                        end do
                else
                        do j=1,dim_of_P(2)
                                if (shocks_under(i)<D(markov(i-1),j)) then
                                        markov(i) = j
                                        exit
                                end if
                        end do
                end if
        end do



        DEALLOCATE (D)
END FUNCTION markov


!Function markov_ind_to_real(shock_space,markov) takes a 1 x T dimensional vector of shocks as indices
!of state space and returns a 1xT vector of shocks translated to the actual realizations
!
!Inputs:
!        - shock_space: 1xM row vector (real)
!        - markov_ind: (1 x T) vector of integers (indices of state space), where T is not fixed.

!Output:
!        - markov_ind_to_real: (1 x T) vector of shock realizations (actual values).
function markov_ind_to_real(shock_space,markov_ind)
        real(kind=dp), dimension(:), intent(in) :: shock_space
        integer, dimension(:) :: markov_ind
        real(kind=dp), dimension(:), allocatable :: markov_ind_to_real

        integer :: T !dimension of vector of shock realizations
        integer :: M !dimension of shock space (number of possible realizations)

        integer :: t_ind !used for cycling

        T = size(markov_ind)
        M = size(shock_space)

        !Allocate memory to output markov_ind_to_real now that we know its dimension
        allocate(markov_ind_to_real(T))

        !First check that the inputs are correct. All indices must be in the range 1,...,M
        if (maxval(markov_ind) > M .or. minval(markov_ind) < 1) then
                stop 'Error: indices in markov_ind out of range.'
        end if

        !Translate the shock realizations from indices to actual realizations
        do t_ind=1,T
                markov_ind_to_real(t_ind) = shock_space(markov_ind(t_ind))
        end do

end function markov_ind_to_real

!compute variance of a series
real(dp) function var_series(ser)
    real(dp), dimension(:), intent(in) :: ser

    real(dp) :: ser_mean !sample mean

    integer :: i

    if(size(ser) < 2) then
        write(*,*) 'Warning: input series must have size at least 2.'
        var_series = 0.0_dp !return variance 0 rather than quit the program altogether
        !to make it easier to locat source of error.
        return
    end if

    ser_mean = 0.0_dp
    do i = 1,size(ser)
        ser_mean = ser_mean + ser(i)/real(size(ser),dp)
    end do

    var_series = 0.0_dp
    do i = 1,size(ser)
        var_series = var_series + ((ser(i) - ser_mean)**2.0_dp)/real(size(ser) - 1,dp)
    end do

end function var_series

!subroutine smallest_indices returns indices of the n smallest elements of a real array x,
!where n is the length of array indices (inout).
subroutine smallest_indices(x,indices)
    real(dp), dimension(:), intent(in) :: x
    integer, dimension(:), intent(inout) :: indices

    real(dp), allocatable, dimension(:) :: x_tmp !copy of x (which we can edit) - so that x can be used with
    !intent (in)

    integer :: min_ind

    !create copy of input vector
    allocate(x_tmp(size(x)))
    x_tmp = x

    do min_ind = 1,size(indices)
        indices(min_ind) = minval(MINLOC(x_tmp)) !minval just to convert integer(1) to integer
        !set the particular value to something very large so that it is
        !not found again as a minimum
        x_tmp(indices(min_ind)) = huge(x) !use intrinsic function for max element?
    end do

    deallocate(x_tmp)

end subroutine smallest_indices

end module mod_tools
