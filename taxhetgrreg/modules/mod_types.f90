module mod_types

use nag_library, only:nag_wp

implicit none

integer, parameter:: sp = selected_real_kind(p=4, r=30)
!integer, parameter:: dp = selected_real_kind(p=8, r=40)

integer, parameter:: dp = nag_wp

integer, parameter:: mp = selected_real_kind(p=6, r=10)

integer, parameter :: long_int_kind = selected_int_kind(r=12)
integer, parameter :: short_int_kind = selected_int_kind(r=6)



end module mod_types
