module mod_work_dist
!Module mod_work_dist contains subroutines/functions which we are used to distribute
!work between images using Coarray Fortran (CAF). This also contains various subroutines
!such as for writing messages, etc.

!General comments
!
!Value function is a multi-dimensional array (4-5 dimensions in my application). To
!distribute work between images (processes), we essentially fold all of these elements
!into a one-dimensional array, then let every image work on a section (stub) of this array.
!When the work is done, we pass all of these stubs to one image where they are combined.
!
!There are many details which will be explained in a separate document! When this is written
!the comments here can be removed.

use mod_types
use mod_taxhetgr

implicit none

!The following type is used for data tranfer of Value function b/w images. We use this because
!to save memory, we want this to be allocated on one image only, and simple
!allocatable coarrays have to be allocated on all images.
!This way, the coarrays will still be allocated on all images
!but on all images except for one all the dimensions used in allocation will be 1
!so pretty much no memory will be used up.
! Now, as long as we didn't
!save any data to the allocated coarrays on the other images we should be fine
!(no memory should be used) but it depends on compiler, perhaps system, etc.
!and it's much safer to do this properly even though the syntax is slightly
!more cumbersome.
!
!also, we could do with just one of V_stub_stacked or V_tmp in principle, but using
!these two allows us to reduce the number of operations needed.
type shared_data_V_type
    real(dp), dimension(:), allocatable :: V_stub_stacked
    real(dp), dimension(:,:,:,:), allocatable :: V_tmp !4 dimensions when I=2
end type shared_data_V_type

!Sharing share of feasible constraints -> define a coarray of reals here, image
!real(sp), codimension[*] :: share_feas_img, then recover these at image 1 and compute weighted average.
!The initial guesses themselves don't need to be accessible to other images.



!Could also use some tricks with coarray of pointers, see the PGAS summer school lecture 2.
contains

!Subroutine get_img_ind serves to compute maximum and minimum index in the folded value
!function that an image will work on. The indices depend on the total number of elements
!of the value function, the index of image that calls the subroutine, and the total number of images.
subroutine get_img_ind(i_img,n_img,V_num_el,min_ind_img,max_ind_img)
    integer, intent(in) :: i_img,n_img,V_num_el
    integer, intent(out) :: min_ind_img, max_ind_img

    integer :: tpi !tasks per image

    if(real(V_num_el)/real(n_img) < 1) then
        write(*,*) 'Error: tasks per image are less than one. Increase the number of grid points or &
        decrease the number of images.'
        error stop
    end if

    tpi = ceiling(real(V_num_el)/real(n_img))
    min_ind_img = 1+(i_img-1)*tpi
    max_ind_img = i_img*tpi
    max_ind_img = min(max_ind_img,V_num_el) !So the image with highest i_img will not exceed the range

end subroutine get_img_ind

!Subroutine send_V_to_spec sends the value function stub from the calling image to the
!image with index img_spec (master image).
subroutine send_V_to_spec(shared_data_V,V_stub,spec_img,min_ind_img,max_ind_img)
    type(shared_data_V_type), codimension[*] :: shared_data_V
    integer, intent(in) :: spec_img,min_ind_img,max_ind_img

    real(dp), dimension(min_ind_img:max_ind_img), intent(in) :: V_stub

    shared_data_V[spec_img]%V_stub_stacked(min_ind_img:max_ind_img) = V_stub(min_ind_img:max_ind_img)
end subroutine send_V_to_spec


!Subroutine get_V_from_spec recovers value function from the special image (spec_img)
!and returns it as V_old (to be saved on the image that called this subroutine)
subroutine get_V_from_spec(shared_data_V,V_old,spec_img)
    type(shared_data_V_type), codimension[*] :: shared_data_V
    integer, intent(in) :: spec_img

    real(dp), dimension(:,:,:,:), intent(out) :: V_old

    V_old(:,:,:,:) = shared_data_V[spec_img]%V_tmp(:,:,:,:)

end subroutine get_V_from_spec


!subroutine get_ind_unfold computes all of the unfolded indices that an image works on, i.e.,
!all of the indices in the multi-dimensional value function which correspond to the indices
!in the vector representation of the value function.
subroutine get_ind_unfold(min_ind_img,max_ind_img,max_ind,ind_unfold_all)
    integer, intent(in) :: min_ind_img,max_ind_img !range of indices in the vector form V that an image works on
    integer, intent(in), dimension(:) :: max_ind !vector of maximum index per state variable

    integer :: i !index for do loop

    integer, dimension(4,min_ind_img:max_ind_img), intent(out) :: ind_unfold_all !4 states (when I=2)

    do i=min_ind_img,max_ind_img
        call ind_FtU(i,max_ind,ind_unfold_all(:,i))
    end do

end subroutine get_ind_unfold

!Subroutine V_unfold recovers value function in the multi-dimensional form
!from the  vector form saved in coarray of shared_data_V (shared_data_V%V_stub_stacked),
!and saves it in shared_data_V%V_tmp. This should only be called on image with
!index spec_img - on all other images, the arrays in the coarray of types will
!be unallocated and it would result in an error.
!
!For some reason, if we directly reshape data from
!the coarray, we get a segmentation error if the value function is even moderately large.
!A solution is to first copy the value function to a temporary array and then reshape it
!and save it.
subroutine V_unfold(shared_data_V)
    type(shared_data_V_type), codimension[*] :: shared_data_V

    real(dp), allocatable, dimension(:) :: V_stub_stacked_loc !Temporary variable so we don't reshape directly

    !Copy data
    allocate(V_stub_stacked_loc(size(shared_data_V%V_stub_stacked)))
    V_stub_stacked_loc = shared_data_V%V_stub_stacked

    !reshape into multi-dimensional representation of V
    shared_data_V%V_tmp = reshape(V_stub_stacked_loc,shape(shared_data_V%V_tmp))

    deallocate(V_stub_stacked_loc)

end subroutine V_unfold

!Subroutine ind_FtU returns indices in the unfolded N-dimensional array corresponding to the
!index ind_folded of the 'folded' array, i.e., of the one-dimensional array which corresponds
!to the multi-dimensional array (using default reshape where order of dimensions is Fortran standard).
!
!This subroutine works for arrays of arbitrary dimension.
subroutine ind_FtU(ind_fold,max_ind,ind_unfold)
    integer, intent(in) :: ind_fold
    integer, intent(in), dimension(:) :: max_ind
    integer, intent(out), dimension(:) :: ind_unfold

    integer :: P,R !Product, Remainder
    integer :: k
    integer :: num_of_dim !number of dimensions of unfolded data

    !Get the product of dimensions
    P=1
    do k=1,size(max_ind)
        P = P*max_ind(k)
    end do

    num_of_dim = size(max_ind)

!Check that the input is correct (can possibly comment this if performance is an issue)
if((num_of_dim /= size(ind_unfold)) .or. any(max_ind<1) .or. ind_fold>P .or. ind_fold < 1) then
    write(*,*) 'Error: wrong input for subroutine ind_FtU!'
    error stop
end if

!If number of dimensions is 1, then ind_unfold is the same as ind_fold then
!stop execution and write an error message - it's not an error as such
!but there is no reason to use this subroutine in this case!
if(num_of_dim == 1) then
    write(*,*) 'num_of_dim in subroutine ind_FtU is 1. No reason to use this, rewrite!'
    error stop
end if

R = ind_fold !'remainder' - i.e., the number of elements that have not been 'placed' yet.
do k=num_of_dim,2,-1
    !when index is k, we recover ind_unfold(k) and prepare
    !remainder to be used in computation of ind_unfold(k-1)
    P = P/max_ind(k) !product of (max_ind(1),...,max_ind(k-1))

    ind_unfold(k) = (R-1)/P + 1

    R = R - (ind_unfold(k)-1)*P
end do
ind_unfold(1) = R !The first index is the remainder. We would get the same result
!if we instead went all the way to k=1, but there would be unnecessary operations.


end subroutine ind_FtU

end module mod_work_dist
