! Copyright (C) 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module model
  use nrtype
  use data, only: nv,ns,nd
  implicit none
  private

  public :: prm,grd

  public :: model_initialize
  public :: model_set_myv
  public :: model_put_myv
  public :: model_collect_prm
  public :: compute_pseudo_likelihood

  public :: mypl
  public :: etot
  public :: fields
  public :: couplings

  ! working variable for the processor
  integer(I4B) :: myv 
  real(DP), allocatable :: fields(:,:) ! NS x NV
  real(DP), allocatable :: couplings(:,:,:) ! NS x NS x NV(NV-1)/2 
  real(DP), allocatable :: p1(:,:) ! NS x NV
  real(DP), allocatable :: p2(:,:,:) ! NS x NS x NV(NV-1)/2 

  ! arrays for fixed residue/sequence
  real(DP), allocatable, save :: my_fields(:) ! NS
  real(DP), allocatable, save :: my_couplings(:,:,:) ! NS x NV x NS
  real(DP), allocatable, save :: my_p1(:) ! NS 
  real(DP), allocatable, save :: my_p2(:,:,:) ! NS x NV x NS
  real(DP), allocatable, save :: my_f1(:) ! NS 
  real(DP), allocatable, save :: my_f2(:,:,:) ! NS x NV x NS

  ! regularization parameters 
  real(DP) :: regularization_strength=0.01_DP ! regularization strength; the default is l2 with regularization_strength=0.01
  integer(I4B) :: regularization_method=2 ! the default is l2; =0 => no regularization 

  ! "cost" functions
  real(DP) :: mypl,ereg,etot

  ! gradients arrays
  real(DP), allocatable, save :: grd(:) ! NS + NS x NV x NS
  real(DP), allocatable, save :: prm(:) ! NS + NS x NV x NS

contains

  subroutine model_initialize(regu,lambda)
    integer(I4B) :: err
    integer(I4B) :: iv,jv
    integer(I4B) :: ind
    integer(I4B), intent(in) :: regu
    real(DP), intent(in) :: lambda


    regularization_method = regu
    regularization_strength = lambda

    allocate(p1(ns,nv),stat=err)
    allocate(p2(ns,ns,nv*(nv-1)/2),stat=err)
    allocate(fields(ns,nv),stat=err)
    allocate(couplings(ns,ns,nv*(nv-1)/2),stat=err)
    allocate(my_f1(ns),stat=err)
    allocate(my_f2(ns,nv,ns),stat=err)
    allocate(my_fields(ns),stat=err)
    allocate(my_couplings(ns,nv,ns),stat=err)
    allocate(my_p1(ns),stat=err)
    allocate(my_p2(ns,nv,ns),stat=err)
    allocate(grd(ns + ns*ns*nv),stat=err)
    allocate(prm(ns + ns*ns*nv),stat=err)


    p1 = 0.0_DP
    p2 = 0.0_DP
    fields = 0.0_DP
    couplings = 0.0_DP       

  end subroutine model_initialize

  subroutine model_zero()
    real(DP), parameter :: small_number=1.e-2_DP
    my_p1 = 0.0_DP
    my_p2 = 0.0_DP
    mypl = 0.0_DP
    etot = 0.0_DP
    select case(regularization_method)
    case(0)
       ereg = 0.0_DP
    case(1)
       write(0,*) "ERROR: l1-regularization not implemented"
       stop
       !       ereg = - regularization_strength * (sum(abs(prm(:ns))) + 0.5_DP * sum(abs(prm(ns+1:))))
       !       ereg = - regularization_strength * small_number * (sum(log(cosh(prm(:ns)/small_number))) + 0.5_DP * sum(log(cosh(prm(ns+1:)/small_number))))
    case(2)
       ereg = - regularization_strength * (sum(prm(:ns)**2) + 0.5_DP * sum(prm(ns+1:)**2))
    case default
       write(0,*) "ERROR: unkown regularization method."
       stop
    end select
  end subroutine model_zero

  subroutine model_set_myv(iv,err) ! my_couplings
    use data, only: nd,data_samples,w
    integer(I4B), intent(in) :: iv
    ! make my_couplings given myv 
    ! must be called before looping on data
    integer(I4B) :: id,jd,jv,k,ind,mys
    integer(I4B) :: err
    integer, allocatable :: list(:)
    real(DP) :: sum,rnd

    myv = iv
    my_p1 = 0.0_DP
    my_p2 = 0.0_DP
    my_f1 = 0.0_DP
    my_f2 = 0.0_DP
    my_fields = 0.0_DP
    my_couplings = 0.0_DP       
    mypl = 0.0_DP
    etot = 0.0_DP
    ereg = 0.0_DP
    prm = 0.0_DP
    grd = 0.0_DP

    ! compute variable-specific arrays of frequencies
    allocate(list(nv),stat=err)
    my_f1 = 0.0_DP
    my_f2 = 0.0_DP
    do id = 1,nd
       list = data_samples(:,id)
       mys = list(myv)
       my_f1(mys) = my_f1(mys) + w(id)
       do jv = 1,nv
          if(jv /= myv) then 
             my_f2(list(jv),jv,mys) = my_f2(list(jv),jv,mys) + w(id)
          end if
       end do
    end do
    deallocate(list)

  end subroutine model_set_myv

  subroutine model_put_myv 
    ! store my_fields in fields 
    ! store my_couplings in couplings
    ! must be called before looping on data
    integer(I4B) :: jv,k,ind

    fields(:,myv) = my_fields
    
    !      --myv--
    !      x x x x
    !   |  1 x x x
    !  jv  2 4 x x
    !   |  3 5 6 x
    !
    ! lower triangle packing: jv > myv

    ! remove gauge before adding to couplings
    call model_gauge()

    fields(:,myv) = my_fields
    do jv = myv+1,nv ! jv > myv
       k = (myv - 1) * nv - myv * (myv + 1) / 2 + jv 
       couplings(:,:,k) = couplings(:,:,k) + my_couplings(:,jv,:)
    end do

    do jv = 1,myv-1 ! myv > jv: submatrices must be transposed 
       k = (jv - 1) * nv - jv * (jv + 1) / 2 + myv
       couplings(:,:,k) =  couplings(:,:,k) + transpose(my_couplings(:,jv,:))
    end do

  end subroutine model_put_myv

  subroutine model_gauge
    integer(I4B) :: jv,is,js
    real(DP) :: mat(ns,ns),arr(ns),marr
    real(DP) :: rsum(ns),csum(ns),totsum

    arr = my_fields
    marr = sum(arr) / real(ns)
    arr = arr - marr
    do jv = 1,nv
       mat = my_couplings(:,jv,:)
       totsum = sum(mat)
       totsum = totsum / real(ns*ns)
       do is = 1,ns
          rsum(is) = sum(mat(is,:))
          csum(is) = sum(mat(:,is))
       end do
       rsum = rsum / real(ns)
       csum = csum / real(ns)
       if(jv /= myv) arr = arr + csum - totsum
       do js = 1,ns
          do is = 1,ns
             mat(is,js) = mat(is,js) - rsum(is) - csum(js) + totsum
          end do
       end do
       my_couplings(:,jv,:) = mat
    end do
    my_fields = arr
    
  end subroutine model_gauge

  subroutine compute_pseudo_conp()
    use data, only: data_samples,w,nd
    integer(I4B) :: list(nv)
    real(DP) :: conp(ns)
    integer(I4B) :: is,jv,kv
    real(DP) :: r,rsum
    real :: start,finish
    real(DP) :: pp,pp0,zz,zz0
    real(DP) :: tmp(nv)
    integer(I4B):: mys
    integer :: ind
    integer :: id
    integer(I4B) :: dd(nv)
    real(DP) :: ww,rnd
    integer(I4B) :: rd

    ! loop over data
    do id = 1,nd
       list = data_samples(:,id)
       ww = w(id)
       mys = list(myv)
       
       ! loop over the states of myv 
       do is = 1,ns
          r = my_fields(is) 
          do jv = 1,nv
             if(myv /= jv) then 
                r = r + my_couplings(list(jv),jv,is) 
             end if
          end do
          conp(is) = exp(r)
       end do
       
       rsum = sum(conp)
       conp = conp / rsum
       
       mypl = mypl + ww * log(conp(mys))
       
       ! update histograms 
       ! loop over the states of myv 
       do is = 1,ns
          pp = conp(is) * ww
          my_p1(is) = my_p1(is) + pp 
          do jv = 1,nv
             if(myv /= jv) then 
                my_p2(list(jv),jv,is) = my_p2(list(jv),jv,is) + pp
             end if
          end do
       end do
       
    end do
    
  end subroutine compute_pseudo_conp

  subroutine compute_pseudo_likelihood(it)
    integer(I4B), intent(in) :: it
    real(DP) :: etot0,de

    call model_parameters_unpack()

    etot0 = etot
    call model_zero()
    call compute_pseudo_conp()
    etot = mypl + ereg
    de = etot-etot0

    call model_parameters_pack()
    
  end subroutine compute_pseudo_likelihood

  subroutine model_parameters_pack
    integer(I4B) :: dim
    real(DP), parameter :: small_number=1.e-2_DP
    integer(I4B) :: is,jv,js

    dim = ns*ns*nv
    ! compute the gradient
    my_p1 = my_p1 - my_f1 
    my_p2 = my_p2 - my_f2 

    select case(regularization_method)
    case(0)
       ! do nothing; no regularization
    case(1) 
       ! l1 regularization
       write(0,*) "ERROR: l1-regularization not implemented"
       stop
       !       my_p1 = my_p1 + regularization_strength * tanh(my_fields/small_number)
       !       my_p2 = my_p2 + 0.5_DP * regularization_strength * tanh(my_couplings/small_number)
    case(2) 
       ! l2 regularization
       my_p1 = my_p1 + 2.0_DP * regularization_strength * my_fields
       my_p2 = my_p2 + 2.0_DP * 0.5_DP * regularization_strength * my_couplings
    case default
       write(0,*) "ERROR: unkown regularization method."
       stop
    end select

    prm(1:ns) = my_fields
    prm(ns+1:) = reshape(my_couplings,(/dim/))

    grd(1:ns) = my_p1
    grd(ns+1:) = reshape(my_p2,(/dim/))
    
  end subroutine model_parameters_pack

  subroutine model_parameters_unpack
    integer :: k

    my_fields = prm(1:ns) 
    my_couplings = reshape(prm(ns+1:),(/ns,nv,ns/))
    
  end subroutine model_parameters_unpack

  subroutine model_collect_prm

    couplings = 0.5_DP * couplings

  end subroutine model_collect_prm
  

end module model
  
