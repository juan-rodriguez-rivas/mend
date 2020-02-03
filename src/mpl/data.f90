! Copyright (C) 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module data
  use nrtype
  implicit none
  private

  public :: nd ! number of data samples
  public :: nv ! number of variables
  public :: ns ! number of states

  public :: data_samples
  public :: w
  public :: neff

  public :: data_read

  integer(I4B) :: nd   ! number of data samples
  integer(I4B) :: nv   ! number of variables
  integer(I4B) :: ns ! number of states per variable

  integer(I4B), allocatable :: data_samples(:,:)
  real(DP), allocatable :: w(:)
  real(DP) :: neff

contains

  subroutine data_initialize(udata)
    use units, only: max_string_length
    use parser, only: parser_nfields
    integer(I4B), intent(in) :: udata
    integer(I4B) :: err
    character(max_string_length) :: line,newline
    integer(I4B) :: nfields
    integer(I4B) :: iv,jv,ind


    ! read the first string
    read(udata,'(a)',iostat=err) line
    call parser_nfields(line,newline,nfields)
    rewind(udata)

    ! set nv
    nv = nfields

    ! count lines for nd
    nd = 0
    do
       read(udata,'(a)',iostat=err) line
       if(err < 0) exit
       nd = nd + 1
    end do
    rewind(udata)

    ! allocate memory for the data

    allocate(data_samples(nv,nd),stat=err)
    allocate(w(nd),stat=err)

    neff = nd

  end subroutine data_initialize

  subroutine data_read(udata,w_id)
    use units, only: max_string_length
    use parser, only: parser_nfields
    integer(I4B), intent(IN) :: udata
    real(DP), intent(IN) :: w_id
    integer(I4B) :: err
    character(max_string_length) :: line,newline
    integer(I4B) :: nfields
    integer(I4B) :: i

    call data_initialize(udata)

    do i = 1,nd
       if(mod(i,100) == 0)write(0,*) 'data: ', i
       read(udata,'(a)',iostat=err) line
       call parser_nfields(line,newline,nfields)
       read(newline,*,iostat=err) data_samples(:,i)
       if(err > 0) write(0,*) 'error: reading data'
    end do

    !    data_samples = data_samples + 1

    if (any(data_samples==0)) then 
       write(0,*) 'ERROR: the data matrix should not contain 0'
       stop
    end if

    ns = maxval(data_samples)

    write(0,*) 'nd: ', nd
    write(0,*) 'nv: ', nv
    write(0,*) 'ns: ', ns
    flush(0)

    if(w_id > 1.E-10_DP) then
       write(0,*) 'computing weights...'
       call data_reweight(w_id)
    else
       w = 1.0_DP / real(nd)
    end if

  end subroutine data_read

  subroutine data_reweight(w_id)
    real(DP), intent(in) :: w_id
    integer(I4B) :: err
    integer(I4B), allocatable :: x(:),y(:)
    integer(I4B) :: id,jd

    allocate(x(nv),y(nv),stat=err)

    w = 1.0_DP
    do id = 1,nd-1
       if(mod(id,1000) == 0) then
          write(0,'(f5.1,a)') 100.0 * real(id) / real(nd),' %'
       end if
       x = data_samples(:,id)
       do jd = id+1,nd
          y = data_samples(:,jd)
          if(count(x == y) >= nint(nv * w_id * 0.01_DP)) then
             w(id) = w(id) + 1.0_DP
             w(jd) = w(jd) + 1.0_DP
          end if
       end do
    end do

    w = 1.0_DP / w
    neff = sum(w)
    w = w / neff
    write(0,*) 'neff: ', neff

  end subroutine data_reweight

end module data
