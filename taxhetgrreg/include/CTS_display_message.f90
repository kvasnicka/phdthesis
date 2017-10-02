!This include file displays message about CTS iteration (number of
!iteration and number of grid points) and saves it in the log file.
if(pars%CTS_split_times>0 .and. i_img == spec_img) then !Only write this message if CTS used
    write(tmp_string1,'(I3)') CTS_ind
    write(tmp_string2,'(I3)') pars%CTS_split_times + 1
    output_string = 'CTS algorithm iteration '//trim(tmp_string1) //' out of '//trim(tmp_string2)
    call double_output(folder_name,trim(output_string))

    write(tmp_string1,'(I3)') N_a
    write(tmp_string2,'(I3)') N_rho
    write(tmp_string3,'(I3)') N_theta

    output_string = 'Number of gridpoints: N_a = '//trim(tmp_string1)//', N_rho = '//trim(tmp_string2)// &
    ', N_theta = '//trim(tmp_string3)
    call double_output(folder_name,trim(output_string))
end if
