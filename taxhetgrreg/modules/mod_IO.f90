module mod_IO
!Module mod_IO contains subroutines used in input/output.

use mod_types

implicit none

contains

subroutine initialize_output(folder_name,parameter_file)
!This subroutine is used to generate name of the folder in which results will be saved.
!Its format is parameter_file_YY_MM_DD_hhmmss. Furthermore, the folder will be created, and the
!parameter file used will be copied into the folder.

use mod_types
 
implicit none

character(256), intent(out) :: folder_name
character*(*), intent(in) :: parameter_file

character(8) :: date
character(10) :: time

character(256) :: shell_cmd

logical :: folder_exists = .false.

integer :: txt_position



call date_and_time(date,time)

!parameter_file already contains .txt suffix and we don't want it to be part of the folder name.
txt_position = index(parameter_file,'.txt')

if (txt_position == 0) then !.txt not present
        folder_name = trim(parameter_file) // '_' // date(3:4) // '_' // date(5:6) //  '_' // date(7:8) // '_' // time(1:6)
else
        folder_name = parameter_file(1:txt_position-1) // '_' // date(3:4) // '_' // date(5:6) &
         //  '_' // date(7:8) // '_' // time(1:6)        
end if

!Test whether the folder exists and stop execution if it does to avoid errors.

!(checking existence of folders problematic with Absoft compiler, so check instead
!whether the log file exist in the folder - we always create it at the same time
!the folder is generated so it's pretty much equivalent unless there was an error in creating the log file).
!With Intel compiler checking existence of folder should work fine - but very low priority change.
inquire(file = 'results/'//trim(folder_name)//'/log.txt',exist = folder_exists)
if (folder_exists) then
        write(*,*) 'Error: folder '//trim(folder_name)//' already exists. Stopping program.'
        error stop
end if


!Also check whether the parameter file exists, if not, stop execution!
inquire(file = 'parameters/'//trim(parameter_file),exist = folder_exists)
if (.not. folder_exists) then
        write(*,*) 'Error: parameter file '//trim(parameter_file)//' does not exist. Stopping execution.'
        error stop
end if

!Create folder:
shell_cmd = 'mkdir -p results/'//trim(folder_name) !-p creates the directory even if directory results doesn't exist
call system(shell_cmd)

!Copy the parameter file into the results folder. This is useful so it's always clear under what
!configuration the results were obtained.
shell_cmd = 'cp parameters/'//trim(parameter_file)//' results/'//trim(folder_name) &
        //'/'//trim(parameter_file)
call system(shell_cmd)


!Create log file
open(unit=21, file = 'results/'//trim(folder_name)//'/log.txt', status = 'replace')
write(21,*) 'Log file (generated '//date(7:8)//'/'//date(5:6)//'/'//date(1:4)//' '&
        //time(1:2)//':'//time(3:4)//':'//time(5:6)//')'
close(unit=21)

end subroutine

!_________________________________________________________________

subroutine double_output(folder_name,output)

!Subroutine double_output saves a string output in file log.txt in folder results\folder_name
!and prints the string output to terminal.
!
!The output is unformatted.

implicit none

character(len=*), intent(in) :: folder_name
character(len=*), intent(in) :: output

logical :: file_exists = .false.

!open file and write the output into it:
inquire(file = 'results/'//trim(folder_name)//'/log.txt',exist = file_exists)
if (.not. file_exists) then
        write(*,*) 'Error: log.txt in folder '//trim(folder_name)//' does not exist. Stopping program.'
        error stop
end if
open(unit=21, file = 'results/'//trim(folder_name)//'/log.txt', status = 'unknown', position = 'append')
write(unit=21, fmt=*)  output                       
close(unit = 21)

!write the same output to console
write(*, fmt=*)  output 

end subroutine



!_________________________________________________________________
!The following function retrieves the name of parameter file from command line
!If there is no name given, the default is chosen (baseline.txt).
character(80) function get_parfilename()
        integer :: narg = 0,  txt_position = 0
        
        narg = command_argument_count()
        
        !default case
        if(narg == 0) then
                get_parfilename = 'baseline.txt'
                return
        end if
        
        if(narg > 1) then
                write(*,*) 'Error: Only one command line argument allowed (name of parameter file).'
                error stop
        end if
        
        !We know that there is only one command line argument - name of parameter file.
        !First check whether the .txt suffix is there. If not, add it. Perverse cases where
        !txt suffix is present more than once or not at the end are neglected.      
        call get_command_argument(1,get_parfilename)
        
        txt_position = index(get_parfilename,'.txt')
        
        if (txt_position == 0) then !.txt not present
                get_parfilename = trim(get_parfilename)//'.txt'
        end if
end function get_parfilename


!_____________________________________________________________________

!The following is terrible (written as one of my first subroutines in Fortran)
!but the priority of rewriting this is very low.

!don't use integers in place of logical vars where
!it attains only 2 values (0,1). Also might get rid of the go to statements.

!Get value of a real parameter (kind = dp) from a text file
!
!Input:
!        - variable_name: name of the variable in the file (character*(*))
!        - file_path: path of the file to read parameters from (character*(*))
!        - parameter_value: the variable which is assigned the value read from the file
!
!Example: GetReal('test','parameters.txt',a) assigns the value 5 to variable a if the file parameters.txt
!contains a row a=5. Comments are also supported in the file (starting with ! or %).
subroutine GetReal(variable_name,file_path,parameter_value)
        character*(*), intent(in) :: variable_name, file_path
        real(kind=dp) :: parameter_value
        
        integer :: file_found
        integer :: ios
        integer :: position
        integer :: variable_found
        integer :: tc1, tc2 !positions of trailing comments
        
        Character(256) :: buffer      
        
        !Assign initial values to variables:
        file_found = 0
        ios = 0
        position = 0
        variable_found = 0
        tc1=0
        tc2=0
        
        !Open the given file
        open(unit=10,action="Read",file=file_path,err=100)
        file_found = 1
        
               
        do while (ios == 0)
                read(10, '(A)', iostat=ios) buffer
                
                !remove leading blanks
                buffer = adjustl(buffer)
                
                !If the first character is % or !, ignore that line (comment)
                if (buffer(1:1) == '%' .or. buffer(1:1) == '!') then    
                        go to 101
                end if
                
                !Find the variable name
                position = index(buffer,variable_name)
                if (position /= 1) go to 101 !If not at the start of the line, go to next line
                
                !find the position of '=' sign (if not there, skip to next line)
                position = index(buffer,'=')
                if (position == 0) go to 101                
                               
                !len_trim(buffer(1:position-1)) is different than length of variable_name, the line
                !does not contain the value of variable which we are looking for (e.g. variable_name='test',
                !and the line contains test1 = 5). Skip the line in this case.
                if(len_trim(buffer(1:position-1)) /= len_trim(variable_name)) go to 101

                !adjustl(buffer(position+1:)) contains everything in the line after '=' with leading blanks shifted to the end
                !Overwrite buffer with this part of the string
                buffer = adjustl(buffer(position+1:))
                
                !Implementation of trailing comments. 
                !Find positions of trailing comments characters ! and %
                tc1 = index(buffer, '!');  tc2 = index(buffer, '%')
                
                if (max(tc1,tc2)==0) then !No trailing comment
                        read(buffer,*, iostat = ios) parameter_value
                else                 
                        if(min(tc1,tc2)==0) then !Only one trailing comment character
                                position = max(tc1,tc2)
                                read(buffer(1:position-1),*, iostat = ios) parameter_value
                        else !Both trailing comment characters found
                                position = min(tc1,tc2)
                                read(buffer(1:position-1),*, iostat = ios) parameter_value                              
                        end if 
                end if
                              
                !If no error (ios == 0), then we found the variable, so increase the counter by 1
                if(ios == 0) variable_found = variable_found + 1
                                               
        101 end do
      
        close(10)
        
!Error handling
100 if(file_found == 0) then
        write(*,*) 'Error: File ', file_path ,' not found.'
        error stop
end if

!If variable_found == 1, everything ok. If it is 0, it was not found, if it is >1, there are conflicts (different values for a variable)
select case(variable_found)
        case(0)
                write(*,*) 'Error: variable ',variable_name, ' not found in file ', file_path, '.'
                error stop
        case(1)
                !everything fine
        case default
                write(*,*) 'Error: variable ',variable_name, ' occurs multiple times in file ', file_path, '.'
                error stop
end select

end subroutine GetReal

!Just like GetReal but for matrices of arbitrary dimension. It is written as a separate function
!so we don't have to declare variables with dimension(1) every time we want to get a value of scalar parameter.
!
!The input parameter_value is now of dimension (:,:). When we want to read a value of row vector, the subroutine
!assumes that these values are on the same row as the name of the vector. For example, the parameter file could
!contain a row: row_vector = 0.2 0.5
!
!If the number of rows of the matrix parameter_value is k>1, the program assumes that these
!values are written on the k rows which follow the name of the variable (no comments are allowed in these k lines).
!For example, the file could contain the following 4 rows:
!maticeA = 
!0.5 0.2 0.3
!0.1 0.2 0.7
!0.5 0.1 0.4


subroutine GetRealMat(variable_name,file_path,parameter_value)
        character*(*), intent(in) :: variable_name, file_path
        real(kind=dp), dimension(:,:) :: parameter_value
        
        integer, dimension(2) :: dim_of_out !dimension of output
        
        integer :: file_found 
        integer :: ios = 0
        integer :: position
        integer :: variable_found
        integer :: tc1, tc2 !positions of trailing comments
        
        Character(256) :: buffer
        
        integer :: i,j !index used in loops
               
        real(kind=dp), dimension(:), allocatable :: row_of_matrix !!used if parameter_value has more than 1 row

        !Assign initial values to variables:
        file_found = 0
        ios = 0
        position = 0
        variable_found = 0
        tc1=0
        tc2=0       
        
        dim_of_out(1) = size(parameter_value,1)
        dim_of_out(2) = size(parameter_value,2)
        
        !Open the given file
        open(unit=10,action="Read",file=file_path,err=102)
        file_found = 1
        
              
        do while (ios == 0)
                read(10, '(A)', iostat=ios) buffer
                
                !remove leading blanks
                buffer = adjustl(buffer)
                
                !If the first character is % or !, ignore that line (comment)
                if (buffer(1:1) == '%' .or. buffer(1:1) == '!') then    
                        go to 105
                end if
                
                !Find the variable name
                position = index(buffer,variable_name)
                if (position /= 1) go to 105 !If not at the start of the line, go to next line
                
                !find the position of '=' sign (if not there, skip to next line)
                position = index(buffer,'=')
                if (position == 0) go to 105                
                               
                !len_trim(buffer(1:position-1)) is different than length of variable_name, the line
                !does not contain the value of variable which we are looking for (e.g. variable_name='test',
                !and the line contains test1 = 5). Skip the line in this case.
                if(len_trim(buffer(1:position-1)) /= len_trim(variable_name)) go to 105
                
                variable_found = variable_found + 1 !If we got here then the variable was found

                if(dim_of_out(1) == 1) then !Handling row vectors

	                !adjustl(buffer(position+1:)) contains everything in the line after '=' with leading blanks shifted to the end
	                !Overwrite buffer with this part of the string
	                buffer = adjustl(buffer(position+1:))
	                
	                !Implementation of trailing comments. 
	                !Find positions of trailing comments characters ! and %
	                tc1 = index(buffer, '!');  tc2 = index(buffer, '%')
	                
	                if (max(tc1,tc2)==0) then !No trailing comment
	                        read(buffer,*, iostat = ios) parameter_value
	                else                 
	                        if(min(tc1,tc2)==0) then !Only one trailing comment character
	                                position = max(tc1,tc2)
	                                read(buffer(1:position-1),*, iostat = ios) parameter_value
	                        else !Both trailing comment characters found
	                                position = min(tc1,tc2)
	                                read(buffer(1:position-1),*, iostat = ios) parameter_value                              
	                        end if 
	                end if
	                              
                
                else !Matrices handled here (more than one row). A matrix with k rows should be given
                !on k lines which follow the row in which the name of the variable is
                 
                !Trailing comments not supported on the lines where matrices are defined in the parameter file (unnecessary complication). Also,
                !the matrix must be given in consecutive rows (no comments in the middle of a matrix allowed). 
                
                allocate(row_of_matrix(dim_of_out(2)))
                do i = 1,dim_of_out(1)
                        read(10, *, iostat=ios) row_of_matrix
                        do j=1,dim_of_out(2)
                                parameter_value(i,j) = row_of_matrix(j)               
                        end do      
                end do
                deallocate(row_of_matrix)

                end if
                                               
        105 end do
      
        close(10)
        
!Error handling
102 if(file_found == 0) then
        write(*,*) 'Error: File ', file_path ,' not found.'
        error stop
end if

!If variable_found == 1, everything ok. If it is 0, it was not found, if it is >1, there are conflicts (different values for a variable)
select case(variable_found)
        case(0)
                write(*,*) 'Error: variable ',variable_name, ' not found in file ', file_path, '.'
                error stop
        case(1)
                !everything fine
        case default
                write(*,*) 'Error: variable ',variable_name, ' occurs multiple times in file ', file_path, '.'
                error stop
        end select                               
end subroutine GetRealMat


!Subroutines GetInt and GetIntMat are exactly the same as GetReal and GetRealMat, the only difference being that
!parameter_value is of integer type. No further description provided.

subroutine GetInt(variable_name,file_path,parameter_value)
        character*(*), intent(in) :: variable_name, file_path
        integer :: parameter_value
        
        integer :: file_found
        integer :: ios
        integer :: position
        integer :: variable_found
        integer :: tc1, tc2 !positions of trailing comments
        
        Character(256) :: buffer      
        
        !Assign initial values to variables:
        file_found = 0
        ios = 0
        position = 0
        variable_found = 0
        tc1=0
        tc2=0
        
        !Open the given file
        open(unit=10,action="Read",file=file_path,err=110)
        file_found = 1
        
               
        do while (ios == 0)
                read(10, '(A)', iostat=ios) buffer
                
                !remove leading blanks
                buffer = adjustl(buffer)
                
                !If the first character is % or !, ignore that line (comment)
                if (buffer(1:1) == '%' .or. buffer(1:1) == '!') then    
                        go to 111
                end if
                
                !Find the variable name
                position = index(buffer,variable_name)
                if (position /= 1) go to 111 !If not at the start of the line, go to next line
                
                !find the position of '=' sign (if not there, skip to next line)
                position = index(buffer,'=')
                if (position == 0) go to 111                
                               
                !len_trim(buffer(1:position-1)) is different than length of variable_name, the line
                !does not contain the value of variable which we are looking for (e.g. variable_name='test',
                !and the line contains test1 = 5). Skip the line in this case.
                if(len_trim(buffer(1:position-1)) /= len_trim(variable_name)) go to 111

                !adjustl(buffer(position+1:)) contains everything in the line after '=' with leading blanks shifted to the end
                !Overwrite buffer with this part of the string
                buffer = adjustl(buffer(position+1:))
                
                !Implementation of trailing comments. 
                !Find positions of trailing comments characters ! and %
                tc1 = index(buffer, '!');  tc2 = index(buffer, '%')
                
                if (max(tc1,tc2)==0) then !No trailing comment
                        read(buffer,*, iostat = ios) parameter_value
                else                 
                        if(min(tc1,tc2)==0) then !Only one trailing comment character
                                position = max(tc1,tc2)
                                read(buffer(1:position-1),*, iostat = ios) parameter_value
                        else !Both trailing comment characters found
                                position = min(tc1,tc2)
                                read(buffer(1:position-1),*, iostat = ios) parameter_value                              
                        end if 
                end if
                              
                !If no error (ios == 0), then we found the variable, so increase the counter by 1
                if(ios == 0) variable_found = variable_found + 1
                                               
        111 end do
      
        close(10)
        
!Error handling
110 if(file_found == 0) then
        write(*,*) 'Error: File ', file_path ,' not found.'
        error stop
end if

!If variable_found == 1, everything ok. If it is 0, it was not found, if it is >1, there are conflicts (different values for a variable)
select case(variable_found)
        case(0)
                write(*,*) 'Error: variable ',variable_name, ' not found in file ', file_path, '.'
                error stop
        case(1)
                !everything fine
        case default
                write(*,*) 'Error: variable ',variable_name, ' occurs multiple times in file ', file_path, '.'
                error stop
end select

end subroutine GetInt


subroutine GetIntMat(variable_name,file_path,parameter_value)
        character*(*), intent(in) :: variable_name, file_path
        integer, dimension(:,:) :: parameter_value
        
        integer, dimension(2) :: dim_of_out !dimension of output
        
        integer :: file_found 
        integer :: ios = 0
        integer :: position
        integer :: variable_found
        integer :: tc1, tc2 !positions of trailing comments
        
        Character(256) :: buffer
        
        integer :: i,j !index used in loops
               
        integer, dimension(:), allocatable :: row_of_matrix !!used if parameter_value has more than 1 row

        !Assign initial values to variables:
        file_found = 0
        ios = 0
        position = 0
        variable_found = 0
        tc1=0
        tc2=0       
        
        dim_of_out(1) = size(parameter_value,1)
        dim_of_out(2) = size(parameter_value,2)
        
        !Open the given file
        open(unit=10,action="Read",file=file_path,err=122)
        file_found = 1
        
              
        do while (ios == 0)
                read(10, '(A)', iostat=ios) buffer
                
                !remove leading blanks
                buffer = adjustl(buffer)
                
                !If the first character is % or !, ignore that line (comment)
                if (buffer(1:1) == '%' .or. buffer(1:1) == '!') then    
                        go to 125
                end if
                
                !Find the variable name
                position = index(buffer,variable_name)
                if (position /= 1) go to 125 !If not at the start of the line, go to next line
                
                !find the position of '=' sign (if not there, skip to next line)
                position = index(buffer,'=')
                if (position == 0) go to 125                
                               
                !len_trim(buffer(1:position-1)) is different than length of variable_name, the line
                !does not contain the value of variable which we are looking for (e.g. variable_name='test',
                !and the line contains test1 = 5). Skip the line in this case.
                if(len_trim(buffer(1:position-1)) /= len_trim(variable_name)) go to 125
                
                variable_found = variable_found + 1 !If we got here then the variable was found

                if(dim_of_out(1) == 1) then !Handling row vectors

	                !adjustl(buffer(position+1:)) contains everything in the line after '=' with leading blanks shifted to the end
	                !Overwrite buffer with this part of the string
	                buffer = adjustl(buffer(position+1:))
	                
	                !Implementation of trailing comments. 
	                !Find positions of trailing comments characters ! and %
	                tc1 = index(buffer, '!');  tc2 = index(buffer, '%')
	                
	                if (max(tc1,tc2)==0) then !No trailing comment
	                        read(buffer,*, iostat = ios) parameter_value
	                else                 
	                        if(min(tc1,tc2)==0) then !Only one trailing comment character
	                                position = max(tc1,tc2)
	                                read(buffer(1:position-1),*, iostat = ios) parameter_value
	                        else !Both trailing comment characters found
	                                position = min(tc1,tc2)
	                                read(buffer(1:position-1),*, iostat = ios) parameter_value                              
	                        end if 
	                end if
	                              
                
                else !Matrices handled here (more than one row). A matrix with k rows should be given
                !on k lines which follow the row in which the name of the variable is
                 
                !Trailing comments not supported on the lines where matrices are defined in the parameter file (unnecessary complication). Also,
                !the matrix must be given in consecutive rows (no comments in the middle of a matrix allowed). 
                
                allocate(row_of_matrix(dim_of_out(2)))
                do i = 1,dim_of_out(1)
                        read(10, *, iostat=ios) row_of_matrix
                        do j=1,dim_of_out(2)
                                parameter_value(i,j) = row_of_matrix(j)               
                        end do      
                end do
                deallocate(row_of_matrix)

                end if
                                               
        125 end do
      
        close(10)
        
!Error handling
122 if(file_found == 0) then
        write(*,*) 'Error: File ', file_path ,' not found.'
        error stop
end if

!If variable_found == 1, everything ok. If it is 0, it was not found, if it is >1, there are conflicts (different values for a variable)
select case(variable_found)
        case(0)
                write(*,*) 'Error: variable ',variable_name, ' not found in file ', file_path, '.'
                error stop
        case(1)
                !everything fine
        case default
                write(*,*) 'Error: variable ',variable_name, ' occurs multiple times in file ', file_path, '.'
                error stop
        end select                               
end subroutine GetIntMat


!subroutine GetString just like GetReal but for strings (character*(*)).
!Don't write ' or " into the parameter file or it will be part of the read string.
subroutine GetString(variable_name,file_path,parameter_value)
        character*(*), intent(in) :: variable_name, file_path
        character(len=*) :: parameter_value
        
        integer :: file_found
        integer :: ios
        integer :: position
        integer :: variable_found
        integer :: tc1, tc2 !positions of trailing comments
        
        Character(256) :: buffer      
        
        !Assign initial values to variables:
        file_found = 0
        ios = 0
        position = 0
        variable_found = 0
        tc1=0
        tc2=0
        
        !Open the given file
        open(unit=10,action="Read",file=file_path,err=100)
        file_found = 1
        
               
        do while (ios == 0)
                read(10, '(A)', iostat=ios) buffer
                
                !remove leading blanks
                buffer = adjustl(buffer)
                
                !If the first character is % or !, ignore that line (comment)
                if (buffer(1:1) == '%' .or. buffer(1:1) == '!') then    
                        go to 101
                end if
                
                !Find the variable name
                position = index(buffer,variable_name)
                if (position /= 1) go to 101 !If not at the start of the line, go to next line
                
                !find the position of '=' sign (if not there, skip to next line)
                position = index(buffer,'=')
                if (position == 0) go to 101                
                               
                !len_trim(buffer(1:position-1)) is different than length of variable_name, the line
                !does not contain the value of variable which we are looking for (e.g. variable_name='test',
                !and the line contains test1 = 5). Skip the line in this case.
                if(len_trim(buffer(1:position-1)) /= len_trim(variable_name)) go to 101

                !adjustl(buffer(position+1:)) contains everything in the line after '=' with leading blanks shifted to the end
                !Overwrite buffer with this part of the string
                buffer = adjustl(buffer(position+1:))
                
                !Implementation of trailing comments. 
                !Find positions of trailing comments characters ! and %
                tc1 = index(buffer, '!');  tc2 = index(buffer, '%')
                
                if (max(tc1,tc2)==0) then !No trailing comment
                        read(buffer,*, iostat = ios) parameter_value
                else                 
                        if(min(tc1,tc2)==0) then !Only one trailing comment character
                                position = max(tc1,tc2)
                                read(buffer(1:position-1),*, iostat = ios) parameter_value
                        else !Both trailing comment characters found
                                position = min(tc1,tc2)
                                read(buffer(1:position-1),*, iostat = ios) parameter_value                              
                        end if 
                end if
                              
                !If no error (ios == 0), then we found the variable, so increase the counter by 1
                if(ios == 0) variable_found = variable_found + 1
                                               
        101 end do
      
        close(10)
        
!Error handling
100 if(file_found == 0) then
        write(*,*) 'Error: File ', file_path ,' not found.'
        error stop
end if

!If variable_found == 1, everything ok. If it is 0, it was not found, if it is >1, there are conflicts (different values for a variable)
select case(variable_found)
        case(0)
                write(*,*) 'Error: variable ',variable_name, ' not found in file ', file_path, '.'
                error stop
        case(1)
                !everything fine
        case default
                write(*,*) 'Error: variable ',variable_name, ' occurs multiple times in file ', file_path, '.'
                error stop
end select

end subroutine GetString
!_________________________________________________________________


!subroutine GetLogical just like GetReal but for logicals.
!Values recognised as .true. are T,1,.true., values recognised as false
!are F,0,.false. Any other value will result in an error.
subroutine GetLogical(variable_name,file_path,parameter_value)
        implicit none

        character*(*), intent(in) :: variable_name, file_path
        logical :: parameter_value

        character(len=256) :: parameter_value_char
        
        integer :: file_found
        integer :: ios
        integer :: position
        integer :: variable_found
        integer :: tc1, tc2 !positions of trailing comments
        
        logical :: found_true = .false., found_false = .false.
        
        Character(256) :: buffer      
        
        !Assign initial values to variables:
        file_found = 0
        ios = 0
        position = 0
        variable_found = 0
        tc1=0
        tc2=0
        
        found_true = .false.; found_false = .false.
        
        !Open the given file
        open(unit=10,action="Read",file=file_path,err=100)
        file_found = 1
        
               
        do while (ios == 0)
                read(10, '(A)', iostat=ios) buffer
                
                !remove leading blanks
                buffer = adjustl(buffer)
                
                !If the first character is % or !, ignore that line (comment)
                if (buffer(1:1) == '%' .or. buffer(1:1) == '!') then    
                        go to 101
                end if
                
                !Find the variable name
                position = index(buffer,variable_name)
                if (position /= 1) go to 101 !If not at the start of the line, go to next line
                
                !find the position of '=' sign (if not there, skip to next line)
                position = index(buffer,'=')
                if (position == 0) go to 101                
                               
                !len_trim(buffer(1:position-1)) is different than length of variable_name, the line
                !does not contain the value of variable which we are looking for (e.g. variable_name='test',
                !and the line contains test1 = 5). Skip the line in this case.
                if(len_trim(buffer(1:position-1)) /= len_trim(variable_name)) go to 101

                !adjustl(buffer(position+1:)) contains everything in the line after '=' with leading blanks shifted to the end
                !Overwrite buffer with this part of the string
                buffer = adjustl(buffer(position+1:))
                
                !Implementation of trailing comments. 
                !Find positions of trailing comments characters ! and %
                tc1 = index(buffer, '!');  tc2 = index(buffer, '%')
                
                if (max(tc1,tc2)==0) then !No trailing comment
                        read(buffer,*, iostat = ios) parameter_value_char
                else                 
                        if(min(tc1,tc2)==0) then !Only one trailing comment character
                                position = max(tc1,tc2)
                                read(buffer(1:position-1),*, iostat = ios) parameter_value_char
                        else !Both trailing comment characters found
                                position = min(tc1,tc2)
                                read(buffer(1:position-1),*, iostat = ios) parameter_value_char                              
                        end if 
                end if
                
                !Now parameter_value_char contains the desired value as a character. We want to convert it to
                !logical value .true. if this is T,.true., or 1, and .false. it this is .false., F, or 0, and
                !display an error message and stop execution of the program otherwise.
                !(already have no leading blanks at this stage and no trailing comments)
                if(index(parameter_value_char,'T')==1 .and. len_trim(parameter_value_char)==1) then
                        found_true = .true.
                end if
                if(index(parameter_value_char,'.true.')==1 .and. len_trim(parameter_value_char)==6) then
                        found_true = .true.
                end if
                if(index(parameter_value_char,'1')==1 .and. len_trim(parameter_value_char)==1) then
                        found_true = .true.
                end if
                
                if(index(parameter_value_char,'F')==1 .and. len_trim(parameter_value_char)==1) then
                        found_false = .true.
                end if
                if(index(parameter_value_char,'.false.')==1 .and. len_trim(parameter_value_char)==7) then
                        found_false = .true.
                end if
                if(index(parameter_value_char,'0')==1 .and. len_trim(parameter_value_char)==1) then
                        found_false = .true.
                end if                
                
                if(found_true) then
                        parameter_value = .true.
                end if
                if(found_false) then
                        parameter_value = .false.
                end if
                
                if(.not. (found_true .or. found_false)) then
                        write(*,*) 'Error: wrong value of variable '//trim(variable_name)//' in parameter file.'
                        write(*,*) '(only values .true.,T,1 or .false.,F,0 allowed.'      
                error stop
                end if
                              
                !If no error (ios == 0), then we found the variable, so increase the counter by 1
                if(ios == 0) variable_found = variable_found + 1
                                               
        101 end do
      
        close(10)
        
!Error handling
100 if(file_found == 0) then
        write(*,*) 'Error: File ', file_path ,' not found.'
        stop
end if

!If variable_found == 1, everything ok. If it is 0, it was not found, if it is >1, there are conflicts (different values for a variable)
select case(variable_found)
        case(0)
                write(*,*) 'Error: variable ',variable_name, ' not found in file ', file_path, '.'
                error stop
        case(1)
                !everything fine
        case default
                write(*,*) 'Error: variable ',variable_name, ' occurs multiple times in file ', file_path, '.'
                error stop
end select

end subroutine GetLogical


end module mod_IO
