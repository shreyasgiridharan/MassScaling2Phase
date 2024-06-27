    module modulefem
    !*********************************************************************
    !    Function: Initialise all variables. Define all allocatables
    !**********************************************************************


    implicit none

    logical, dimension(:), allocatable :: IsFixDof
    integer :: ndtn = 0, ITime = 0
    integer :: NDivX = 0, NDivY = 0
    integer :: NEl = 0, NNod = 0
    integer, dimension(:,:), allocatable :: ICon
    double precision :: gx = 0.D0, gy = 0.D0, rhoD, rhoS, rhoSat, poro, rhoW=1.d0, kappa, PWaveS
    double precision :: Meshdx = 0.D0, Meshdy = 0.D0
    double precision, dimension(:), allocatable :: GrvF, InrF, MasS, Area, Areai, ExtF, v0, V, dis, MasW, MasbW, InrFW, ExtFW, DragF, GrvFW, ExtP, W, w0, InrF0, InrFW0, RhoSol, RhoFlu, Perm, Damp, DampW
    double precision, dimension(:,:), allocatable :: NodCo
    double precision, dimension(:,:,:), allocatable :: SigG, SigE, EpsG, F, HS, PoreP, TempP,  Plastind, FBARNEW, dEpsG
    double precision, dimension(:,:,:,:), allocatable :: B
    double precision :: MatProp(7), dt = 0.D0, FBAR, epsv, Props(20)
    double precision :: BulkW = 2200000.d0 ! Bulk Modulus of Water 2200000 kN/m2
    double precision :: PWave = 0.d0 !
    double precision :: MinMeshSize = 0.d0 !
    double precision :: dtime = 1.d0 !
    double precision :: dtcr= 1.05d0 !
    double precision :: time = 0.d0
    double precision :: totaltime = 0.d0
    double precision :: Beta_V = 100.d0
    double precision :: alphaS=0.05d0
    double precision :: alphaW=0.05d0

    integer :: nplast = 0.d0
    double precision :: PHI, PSI, COH

    double precision :: m_n160, m_kge, m_kb, m_me, m_ne, &
        m_kgp, m_np, m_phipt, m_phif, &
        m_hfac1, m_hfac2, m_hfac3, m_hfac4, m_hfac5, m_hfac6, &
        m_rf, m_anisofac, m_coh, pa, m_ten, output

    integer, parameter :: DATUnit = 2
    integer, parameter :: LOGUnit = 2
    integer, parameter :: MSHUnit = 3
    integer, parameter :: RESUnit = 4

    contains

    subroutine solve()
    !**********************************************************************
    !    Function: Compute solution over a loop
    !**********************************************************************

    implicit none
    integer :: IIter, IStep = 0, IPTime, dPsi
    integer :: iprint, imeshwrite, iter, niter, Id, I, IEl, Ig, INod(4)
    double precision :: duration, factor, error1, error2, tolerance, sum1, sum2, sum3, sum4, error, error0
    real :: start, finish
    iprint = 1000
    imeshwrite = 1000
    factor=0.d0
    niter = 100
    tolerance=0.0001d0
    error1=9999.d0; error2=9999.d0
    error=9999.d0
    open(LOGUnit, file = 'output.log')
    call initial()
    call WtMSH(MSHUnit)
    call CPU_TIME(start)

    duration  = totaltime/ndtn
    time=0.d0

    ! Run the solution over a loop
    do ITime = 1, ndtn ! physical time scale
        time = time + duration
        factor=ITime/ndtn
        iter = 0

        RhoFlu=0.d0; RhoSol=0.d0
        call DetRho()
        call Map2nod()

        ! Compute the error in the substep and iterate till convergance
        ! Prof. Stolle is not really sure on which approach is more efficient
        ! Will have to perform some studies, BVL to evaluate
        ! Here, implementing both types
        do while(iter.le.niter.and.error.gt.tolerance)
            !do while(error.gt.tolerance)
            !do while(error2.gt.tolerance)
            iter = iter + 1
            !call Map2NodIntr()
            inrf0=inrf; inrfw0=inrfw
            call update_small()
            call Map2NodIntr()
            sum1=0.d0; sum2=0.d0; sum3=0.d0; sum4=0.d0;
            !do I = 1, NNod
            !    Id = (I - 1) * 2
            !    if (.not.IsFixDof(Id + 1).and.(MasS(Id + 1) .gt. 0.D0)) then
            !        sum1 = sum1 + (GrvF(Id+1)-InrF(Id+1))**2.d0 + (GrvFW(Id+1)-InrFW(Id+1))**2.d0
            !        sum2 = sum2 + (GrvF(Id+1))**2.d0 + (GrvFW(Id+1))**2.d0
            !    end if
            !    if (.not.IsFixDof(Id + 2).and.(MasS(Id + 2) .gt. 0.D0)) then
            !        sum1 = sum1 + (GrvF(Id+2)-Inrf(Id+2))**2.d0 + (GrvFW(Id+2)-InrFW(Id+2))**2.d0
            !        sum2 = sum2 + (GrvF(Id+2))**2.d0 + (GrvFW(Id+2))**2.d0
            !    end if
            !end do

            !sum1 = norm2(Damp) + norm2(DampW)
            !sum2 = norm2(GrvF) + norm2(GrvFW)

            sum1 = norm2(inrf+inrfw - inrf0-inrfw0)
            sum2 = norm2(inrf+inrfw)

            !do I = 1, NNod
            !    Id = (I - 1) * 2
            !    if (.not.IsFixDof(Id + 1).and.(MasW(Id + 1) .gt. 0.D0)) then
            !        sum3 = sum3 + (GrvFW(Id+1)-InrFW(Id+1))**2.d0
            !        sum4 = sum4 + (GrvFW(Id+1))**2.d0
            !    end if
            !    if (.not.IsFixDof(Id + 2).and.(MasW(Id + 2) .gt. 0.D0)) then
            !        sum3 = sum3 + (GrvFW(Id+2)-InrFW(Id+2))**2.d0
            !        sum4 = sum4 + (GrvFW(Id+2))**2.d0
            !    end if
            !end do

            !error1=sqrt(sum1/sum2)!*ndtn
            !error2=sqrt(sum3/sum4)!*ndtn
            !error=max(error1,error2)

            error=sqrt(sum1/sum2)!*ndtn

            !tempp=0.d0
            !do IEl = 1, NEl
            !    INod(:) = ICon(:, IEl)
            !    do ig = 1, 4
            !        do I = 1, 4
            !            Id = (INod(I) - 1) * 2
            !            If(MasS(ID+1).gt.0.d0) then
            !                tempp(1:3,ig,iel) = porep(1:3,ig,iel)*MasS(Id +1)
            !            end if
            !            If(MasS(ID+1).gt.0.d0) then
            !                PoreP(1:3,ig,iel) = tempp(1:3,ig,iel)/MasS(Id+1)
            !            end if
            !        end do
            !        ! Strain increment
            !    enddo !gauss
            !end do


        end do

        !if(error0<error) then
        !    niter = niter + 0.2*niter
        !else
        !    niter = niter - 0.1*niter
        !end if

        !!if(niter>10000) niter = 10000
        !
        error0=error

        !dis = dis + v * dtime
        !v0 = v
        !w0 = w

        !v0=w0

        time = time + totaltime

        if(itime.eq.1 .or. (itime/imeshwrite)*imeshwrite == itime .or. itime == ndtn) then
            IStep = IStep + 1
            call WtRES(IStep)
        endif
        if(itime.eq.1 .or. (itime/iprint)*iprint == itime .or. itime == ndtn) then
            write (*, *) sigg(2,1,1), itime*dt
            write (LOGUnit, *) PoreP(2,1,1), error, iter
        endif
    enddo
    call CPU_TIME(finish)
    write(*,*) 'Time = ',finish-start,' seconds'
    write(*,*) 'Time = ',(finish-start)/60.d0,' minutes'
    write(*,*) 'Time = ',(finish-start)/3600.d0,' hours'
    !read(*,*)
    close(2)

    end subroutine solve

    subroutine readdata()
    !**********************************************************************
    !    Function: Read the file
    !**********************************************************************

    implicit none
    double precision :: r1(2), r2(2), LPos(2), temp, Gpy, PI, minx,miny
    integer :: I, J, IEl, INod(4), ig

    open(DATUnit, file = 'input.dat')
    read(DATUnit, *) NDivX, NDivY
    read(DATUnit, *) Meshdx, Meshdy
    NNod = (NDivX + 1)*(NDivY + 1)
    NEl = NDivX * NDivY
    allocate (NodCo(2, NNod)); nodco = 0.d0
    allocate (ICon(4, NEl)); icon = -1
    allocate (IsFixDof(2 * NNod))
    allocate(MasS(2 * NNod)); MasS = 0.d0
    allocate(GrvF(2 * NNod)); GrvF = 0.d0
    allocate(InrF(2 * NNod)); InrF = 0.d0
    allocate(InrF0(2 * NNod)); InrF0 = 0.d0
    allocate(ExtF(2 * NNod)); ExtF = 0.d0

    allocate(ExtP(2 * NNod)); ExtP = 0.d0
    allocate(InrFW(2 * NNod)); InrFW = 0.d0
    allocate(InrFW0(2 * NNod)); InrFW0 = 0.d0
    allocate(GrvFW(2 * NNod)); GrvFW = 0.d0
    allocate(MasW(2 * NNod)); MasW = 0.d0
    allocate(MasbW(2 * NNod)); MasbW = 0.d0
    allocate(DragF(2 * NNod)); DragF = 0.d0

    allocate(Damp(2 * NNod)); Damp = 0.d0
    allocate(DampW(2 * NNod)); DampW = 0.d0


    allocate(v(2 * NNod)); V = 0.d0
    allocate(v0(2 * NNod)); V0 = 0.d0
    allocate(W(2 * NNod)); W = 0.d0
    allocate(W0(2 * NNod)); W0 = 0.d0

    allocate(dis(2 * NNod)); dis = 0.d0
    allocate(B(2, 4, NEl, 4)); B = 0.d0
    allocate(Area(NEl)); Area = 0.d0
    allocate(Areai(NEl)); Areai = 0.d0
    allocate(Sigg(3, 4, NEl)); Sigg = 0.d0
    allocate(SigE(3, 4, NEl)); SigE = 0.d0
    allocate(EpsG(3, 4, NEl)); EpsG = 0.d0
    allocate(dEpsG(3, 4, NEl)); dEpsG = 0.d0
    allocate(PoreP(3, 4, NEl)); PoreP = 0.d0

    allocate(TempP(3, 4, NEl)); TempP = 0.d0

    allocate(RhoSol(NEl)); RhoSol = 0.d0
    allocate(RhoFlu(NEl)); RhoFlu = 0.d0
    allocate(Perm(NEl)); Perm = 0.d0

    allocate(HS(4, NEL, 4)); HS = 0.d0

    allocate(PlastInd(1,4,Nel)); PlastInd = 0.d0 ! State variable
    allocate(FBARNEW(1, 4, NEl));  FBARNEW = 0.d0


    do I = 1, NNod
        read(DATUnit, *) temp, NodCo(1, I), NodCo(2, I), temp
    end do

    do I = 1, NEl
        read(DATUnit, *) temp, ICon(1, I), ICon(2, I), ICon(3, I), ICon(4, I)
    end do

    read(DATUnit, *) MatProp(1), MatProp(2), MatProp(3), MatProp(4), MatProp(5) !E,nu,rhoS,n,kappa(conductivity)
    read(DATUnit, *) gx, gy
    read(DATUnit, *) ndtn, dt
    close(1)
    Bulkw = 2.2E6
    PWave = MatProp(1)*(1.d0-MatProp(2))/((1.d0+MatProp(2))*(1.d0-2.d0*MatProp(2)))

    do IEl = 1, NEl
        INod(:) = ICon(:, IEl)
        Area(IEl) = abs(NodCo(1,INod(3))-NodCo(1,INod(1)))* abs(NodCo(2,INod(3))-NodCo(2,INod(1)))
    end do

    MinMeshSize=999999.d0; minx=99999.d0; miny=99999.d0
    do I=1,NNod-1
        MinX=abs(NodCo(1,I)-NodCo(1,I+1))
        MinY=abs(NodCo(2,I)-NodCo(2,I+1))
        if ((MinMeshSize.gt.MinX).and.(MinX.gt.0.0001d0)) MinMeshSize=MinX
        if ((MinMeshSize.gt.MinY).and.(MinY.gt.0.0001d0)) MinMeshSize=MinY
    end do

    do I=1,nnod
        NodCo(1,I) = NodCo(1,I) * Meshdx
        NodCo(2,I) = NodCo(2,I) * Meshdy
    enddo
    PI = 4.d0 * atan(1.d0) ! Approximation for PI
    !PHI  = SIN(PI*MatProp(5)/180.D0)
    !PSI  = SIN(PI*MatProp(6)/180.D0)
    !COH = MatProp(7)*COS(PI*PHI/180.D0)

    ! Apply boundary conditions
    IsFixDof = .false.
    !IsFixDof(1) = .true.
    !IsFixDof(2) = .true.
    !IsFixDof(3) = .true.
    !IsFixDof(4) = .true.
    !IsFixDof(5) = .true.
    !IsFixDof(7) = .true.

    do I = 1, NNod ! vertical bar problem
        if ((NodCo(2, I) .eq. 0.d0)) then ! bottom
            IsFixDof((I - 1) * 2 + 2) = .true.
            IsFixDof((I - 1) * 2 + 1) = .true.
        end if
        if ((NodCo(1, I) .eq. 0.d0)) then ! left
            IsFixDof((I - 1) * 2 + 1) = .true.
            !IsFixDof((I - 1) * 2 + 2) = .true.
        end if
        if ((NodCo(1, I) .eq. 1.d0)) then ! right
            IsFixDof((I - 1) * 2 + 1) = .true.
            !IsFixDof((I - 1) * 2 + 2) = .true.
        end if
    end do

    !Initial stresses, gravity stress initialisation
    !do IEl = 1, nel
    !    INod(:) = ICon(:, IEl)
    !    !do ig = 1, 4
    !        GPy=(NodCo(2,INod(1))+NodCo(2,INod(4)))/2.d0
    !        SigE(2,1:4,IEl)=(0.5d0-GPy)*((1.d0-MatProp(4))*MatProp(3)+MatProp(4))*gy
    !        SigE(1,1:4,IEl)=SigE(2,1:4,IEl)/2.0
    !        !SigE(3,1:4,IEl)=(SigE(1,1:4,IEl)+SigE(2,1:4,IEl))/2.0
    !        PoreP(1:2,1:4,IEl) = (0.5d0-GPy)*gy
    !    !enddo !gauss
    !end do !nel
    !SigG=SigE+PoreP


    end subroutine readdata

    subroutine Initial()

    !**********************************************************************
    !    Function: Calculates B-Strain Displacement matrix (1 and 4 gauss points)
    !**********************************************************************

    implicit none

    integer :: IEl, I, J, K, INod(4), Id, ig
    double precision :: LPos(2, 1), dNxi(2, 4), Ja(2, 2), JaI(2, 2), A
    double precision :: xi, eta, rp, rm, sp, sm, temp

    !Porosity and Soild density
    rhoS = MatProp(3)
    Poro = MatProp(4)
    kappa = MatProp(5) ! conductivity or Darcy Permeability

    !dry and saturated density
    rhoD= (1.d0-poro)*rhoS
    rhoSat= rhoD + poro * rhoW

    temp = 1.d0/sqrt(3.d0) ! 4 Gauß Points
    !temp=0.d0 ! 1 Gauß Points

    B = 0.0d0

    do IEl = 1, nel
        INod(:) = ICon(:, IEl)

        do ig = 1, 4

            select case (ig)
            case(1)
                xi = -temp
                eta = -temp

            case (2)
                xi = temp
                eta = -temp

            case(3)
                xi = temp
                eta = temp

            case(4)
                xi = -temp
                eta = temp
            end select

            rp = 1.0 + xi
            rm = 1.0 - xi;
            sp = 1.0 + eta;
            sm = 1.0 - eta;

            dNxi(1, 1) = -0.25D0 * sm; dNxi(1, 2) = +0.25D0 * sm; dNxi(1, 3) = +0.25D0 * sp
            dNxi(1, 4) = -0.25D0 * sp; dNxi(2, 1) = -0.25D0 * rm; dNxi(2, 2) = -0.25D0 * rp
            dNxi(2, 3) = +0.25D0 * rp; dNxi(2, 4) = +0.25D0 * rm

            HS(1, iel, ig) = (1.D0 - xi)*(1.D0 - eta)/4.D0 ! Shape Functions
            HS(2, iel, ig) = (1.D0 + xi)*(1.D0 - eta)/4.D0
            HS(3, iel, ig) = (1.D0 + xi)*(1.D0 + eta)/4.D0
            HS(4, iel, ig) = (1.D0 - xi)*(1.D0 + eta)/4.D0

            Area(IEl) = abs(NodCo(1, INod(3)) - NodCo(1, INod(1))) * abs(NodCo(2, INod(3)) - NodCo(2, INod(1)))

            Ja = 0.0D0

            do I = 1, 2
                do J = 1, 2
                    do K = 1, 4
                        Ja(I, J) = Ja(I, J) + dNxi(I, K) * NodCo(J, INod(K)) ! Jacobian
                    end do
                end do
            end do
            A = Ja(1, 1) * Ja(2, 2) - Ja(1, 2) * Ja(2, 1)

            if (A .gt. 0.D0) then
                JaI(1, 1) = +Ja(2, 2)/A; JaI(1, 2) = -Ja(1, 2)/A ! Inverse of Jacobian
                JaI(2, 1) = -Ja(2, 1)/A; JaI(2, 2) = +Ja(1, 1)/A
            else
                write(LOGUnit, *) 'negative or zero Jacobian !!'; stop
            end if

            do J = 1, 4
                B(1, J, IEl, ig) = dNxi(1, J) * JaI(1, 1) + dNxi(2, J) * JaI(1, 2) ! Strain displacement matrix
                B(2, J, IEl, ig) = dNxi(1, J) * JaI(2, 1) + dNxi(2, J) * JaI(2, 2)
            end do
        enddo
    enddo

    ! call Map2Nod()

    end subroutine Initial


    subroutine DetRho()
    !**********************************************************************
    !    Function: Determine critical densities. Refer to Prof. Stolle's notes on Dynamic Relaxation
    !**********************************************************************

    implicit none
    integer :: IEl, INod(4)
    double precision :: cd, RhTilda, RhSat, RhSolid, RhTildaW, RhSatW, RhSolidW, Rh, RhFluid, RhPerm

    cd = 0.d0; RhoSol=0.d0; RhTilda = 0.d0; RhSat=0.d0; RhSolid=0.d0; RhTildaW = 0.d0; RhSatW=0.d0; RhSolidW=0.d0; Rh=0.d0; RhPerm = 0.d0; Perm = 0.d0; RhFluid=0.d0
    do IEl = 1, NEl
        INod(:) = ICon(:, IEl)
        cd = MinMeshSize / dtcr
        PWaveS = MatProp(1)*(1.d0-MatProp(2))/((1.d0+MatProp(2))*(1.d0-2.d0*MatProp(2)))
        RhSat = (PWaveS + BulkW/poro)/(cd**2.d0)
        RhTilda = PWaveS/(cd**2.d0)
        !RhFluid = (RhTilda - RhSat) / (1.d0/poro - 2.d0)
        RhPerm = dtcr*RhFluid*abs(gy)/(2.d0*RhTilda)
        !RhFluid = (BulkW/poro)/(cd**2.d0)
        RhFluid = (BulkW)/(cd**2.d0)

        !RhSolid = (RhSat - MatProp(4)*rhoW) / (1.d0 - MatProp(4))

        !RhFluid = RhSat * (poro*BulkW/(poro*PWaveS + BulkW))
        RhSolid = (RhSat - MatProp(4)*RhFluid) / (1.d0 - MatProp(4))


        RhoSol(IEl) = RhoSol(IEl) + RhSat*Area(iel)/4.d0
        RhoFlu(IEl) = RhoFlu(IEl) + RhFluid*Area(iel)/4.d0
        Perm(IEl) = Perm(IEl) + RhPerm*Area(iel)/4.d0
    end do

    do IEl = 1, NEl
        RhoSol(IEl) = RhoSol(IEl)/Area(iel)*4.d0
        RhoFlu(IEl) = RhoFlu(IEl)/Area(iel)*4.d0
        Perm(IEl) = Perm(IEl)/Area(iel)*4.d0
    end do

    end subroutine DetRho


    subroutine Map2Nod()

    !**********************************************************************
    !    Function: Internal and External Forces, Mass at each nodes
    !**********************************************************************

    implicit none
    integer :: I, IEl, Id, INod(4), ig, J, count
    double precision :: factor
    GrvF = 0.D0; InrF = 0.D0; ExtF = 0.D0; MasS = 0.D0
    GrvFW = 0.D0; InrFW = 0.D0; ExtP = 0.D0; MasW = 0.D0; MasbW = 0.d0; DragF = 0.d0; DampW = 0.d0; Damp = 0.d0

    !if(itime*dt.gt.5000.d0) then
    !    factor = 1.d0
    !else
    !    factor = itime * dt / 5000.d0
    !endif
    !write(*,*) itime*dt,factor

    if(itime*dt.gt.5000.d0) then; alphaS=0.15d0;alphaW=0.15d0; end if


        do IEl = 1, nel
            INod(:) = ICon(:, IEl)
            do ig = 1, 4 ! gauss
                do I = 1, 4
                    Id = (INod(I) - 1) * 2
                    ! Soil Part
                    InrF(Id + 1) = InrF(Id + 1)+(Sigg(1,ig,IEl) * B(1,I,iel,ig) + Sigg(3,ig,IEl) * B(2, I, iel, ig)) * Area(iel)/4.d0
                    InrF(Id + 2) = InrF(Id + 2)+(Sigg(3,ig,IEl) * B(1,I,iel,ig) + Sigg(2,ig,IEl) * B(2, I, iel, ig)) * Area(iel)/4.d0
                    GrvF(Id + 1) = GrvF(Id + 1) + Area(iel) * rhoSat * hs(i, iel, ig) * gx /4.d0 !* factor
                    GrvF(Id + 2) = GrvF(Id + 2) + Area(iel) * rhoSat * hs(i, iel, ig) * gy /4.d0 !* factor
                    MasS(Id + 1) = MasS(Id + 1) + Area(iel) * (1 - poro) * rhoSol(IEl) * hs(i, iel, ig) /4.d0
                    MasS(Id + 2) = MasS(Id + 2) + Area(iel) * (1 - poro) * rhoSol(IEl) * hs(i, iel, ig) /4.d0
                    ! Water part
                    InrFW(Id + 1) = InrFW(Id + 1) + (PoreP(1,ig,IEl) * B(1,I,iel,ig)) * Area(iel)/4.d0
                    InrFW(Id + 2) = InrFW(Id + 2) + (PoreP(2,ig,IEl) * B(2,I,iel,ig)) * Area(iel)/4.d0
                    GrvFW(Id + 1) = GrvFW(Id + 1) + Area(iel) * rhoW * hs(i, iel, ig) * gx /4.d0 !* factor
                    GrvFW(Id + 2) = GrvFW(Id + 2) + Area(iel) * rhoW * hs(i, iel, ig) * gy /4.d0 !* factor
                    MasW(Id + 1)  = MasW(Id + 1) + Area(iel) *  RhoFlu(IEl) * hs(i, iel, ig) /4.d0
                    MasW(Id + 2)  = MasW(Id + 2) + Area(iel) *  RhoFlu(IEl) * hs(i, iel, ig) /4.d0
                    MasbW(Id + 1) = MasbW(Id + 1) + Area(iel) * (poro) * RhoFlu(IEl) * hs(i, iel, ig) /4.d0
                    MasbW(Id + 2) = MasbW(Id + 2) + Area(iel) * (poro) * RhoFlu(IEl) * hs(i, iel, ig) /4.d0
                    DragF(Id + 1) = DragF(Id + 1) + ( Area(iel) *  rhoW * poro * 9.81d0 * hs(i, iel, ig) / (kappa * 4.d0)) * (W0(Id+1) - V0(Id+1))
                    DragF(Id + 2) = DragF(Id + 2) + ( Area(iel) *  rhoW * poro * 9.81d0 * hs(i, iel, ig) / (kappa * 4.d0)) * (W0(Id+2) - V0(Id+2))
                end do
            enddo !gauss
        end do !nel

        do IEl = 1, nel
            INod(:) = ICon(:, IEl)
            do ig = 1, 4 ! gauss
                do I = 1, 4
                    Id = (INod(I) - 1) * 2
                    if(abs(W0(Id + 1)).gt.0.d0) DampW(Id + 1) = DampW(Id + 1) - alphaW * abs(GrvFW(Id + 1) - InrFW(Id + 1))*sign(1.d0,W0(Id + 1))
                    if(abs(W0(Id + 2)).gt.0.d0) DampW(Id + 2) = DampW(Id + 2) - alphaW * abs(GrvFW(Id + 2) - InrFW(Id + 2))*sign(1.d0,W0(Id + 2))
                end do
            enddo !gauss
        end do !nel

        do IEl = 1, nel
            INod(:) = ICon(:, IEl)
            do ig = 1, 4 ! gauss
                do I = 1, 4
                    Id = (INod(I) - 1) * 2
                    if(abs(V0(Id + 1)).gt.0.d0) Damp(Id + 1) = Damp(Id + 1) - alphaS * abs((GrvF(Id + 1) - InrF(Id + 1)) - (GrvFW(Id + 1) - InrFW(Id + 1)))*sign(1.d0,V0(Id + 1))
                    if(abs(V0(Id + 2)).gt.0.d0) Damp(Id + 2) = Damp(Id + 2) - alphaS * abs((GrvF(Id + 2) - InrF(Id + 2)) - (GrvFW(Id + 2) - InrFW(Id + 2)))*sign(1.d0,V0(Id + 2))
                end do
            enddo !gauss
        end do !nel


    end subroutine Map2Nod

    subroutine Map2NodIntr()

    !**********************************************************************
    !    Function: Internal forces  at each nodes
    !**********************************************************************

    implicit none
    integer :: I, IEl, Id, INod(4), ig, J
    double precision :: factor
    InrF = 0.D0; InrFW=0.d0; Damp=0.d0; DampW=0.d0; !DragF=0.d0

    if(itime*dt.gt.5000.d0) then; alphaS=0.15d0;alphaW=0.15d0; end if


        do IEl = 1, nel
            INod(:) = ICon(:, IEl)
            do ig = 1, 4 ! gauss
                do I = 1, 4
                    Id = (INod(I) - 1) * 2
                    InrF(Id + 1) = InrF(Id + 1) + (Sigg(1,ig,IEl) * B(1,I,iel,ig) + Sigg(3,ig,IEl) * B(2, I, iel, ig)) * Area(iel)/4.d0
                    InrF(Id + 2) = InrF(Id + 2) + (Sigg(3,ig,IEl) * B(1,I,iel,ig) + Sigg(2,ig,IEl) * B(2, I, iel, ig)) * Area(iel)/4.d0
                    InrFW(Id + 1) = InrFW(Id + 1) + (PoreP(1,ig,IEl) * B(1,I,iel,ig)) * Area(iel)/4.d0
                    InrFW(Id + 2) = InrFW(Id + 2) + (PoreP(2,ig,IEl) * B(2,I,iel,ig)) * Area(iel)/4.d0

                    !DragF(Id + 1) = DragF(Id + 1) + ( Area(iel) *  rhoW * poro * 9.81d0 * hs(i, iel, ig) / (kappa * 4.d0)) * (W0(Id+1) - V0(Id+1))
                    !DragF(Id + 2) = DragF(Id + 2) + ( Area(iel) *  rhoW * poro * 9.81d0 * hs(i, iel, ig) / (kappa * 4.d0)) * (W0(Id+2) - V0(Id+2))

                    !DampW(Id + 1) = DampW(Id + 1) - alphaW * abs(GrvFW(Id + 1) - InrFW(Id + 1))*sign(1.d0,W(Id + 1))
                    !DampW(Id + 2) = DampW(Id + 2) - alphaW * abs(GrvFW(Id + 2) - InrFW(Id + 2))*sign(1.d0,W(Id + 2))

                    !Damp(Id + 1) = Damp(Id + 1) - alphaS * (abs(GrvF(Id + 1) - InrF(Id + 1)) - abs(GrvFW(Id + 1) - InrFW(Id + 1)))*sign(1.d0,V(Id + 1))  + DampW(Id + 1)
                    !Damp(Id + 2) = Damp(Id + 2) - alphaS * (abs(GrvF(Id + 2) - InrF(Id + 2)) - abs(GrvFW(Id + 2) - InrFW(Id + 2)))*sign(1.d0,V(Id + 2))  + DampW(Id + 2)

                end do
            enddo !gauss
        end do !nel


        do IEl = 1, nel
            INod(:) = ICon(:, IEl)
            do ig = 1, 4 ! gauss
                do I = 1, 4
                    Id = (INod(I) - 1) * 2
                    if(abs(W0(Id + 1)).gt.0.d0) DampW(Id + 1) = DampW(Id + 1) - alphaW * abs(GrvFW(Id + 1) - InrFW(Id + 1))*sign(1.d0,W0(Id + 1))
                    if(abs(W0(Id + 2)).gt.0.d0) DampW(Id + 2) = DampW(Id + 2) - alphaW * abs(GrvFW(Id + 2) - InrFW(Id + 2))*sign(1.d0,W0(Id + 2))
                end do
            enddo !gauss
        end do !nel

        do IEl = 1, nel
            INod(:) = ICon(:, IEl)
            do ig = 1, 4 ! gauss
                do I = 1, 4
                    Id = (INod(I) - 1) * 2
                    if(abs(V0(Id + 1)).gt.0.d0) Damp(Id + 1) = Damp(Id + 1) - alphaS * abs((GrvF(Id + 1) - InrF(Id + 1)) - (GrvFW(Id + 1) - InrFW(Id + 1)))*sign(1.d0,V0(Id + 1))
                    if(abs(V0(Id + 2)).gt.0.d0) Damp(Id + 2) = Damp(Id + 2) - alphaS * abs((GrvF(Id + 2) - InrF(Id + 2)) - (GrvFW(Id + 2) - InrFW(Id + 2)))*sign(1.d0,V0(Id + 2))
                end do
            enddo !gauss
        end do !nel




    end subroutine Map2NodIntr


    subroutine  Update_small()
    !**********************************************************************
    !    Function: Small deformation formulation
    !**********************************************************************

    implicit none
    integer :: IEl, INod(4), Id, I, ig, J
    double precision :: delV(4), deps(3), Bulk = 2e6, delW(4), dpore, aW, temp=0.d0, dampf=0.001d0, dSig(4)
    !Initial Velocity
    V = 0.d0
    W = 0.d0
    do I = 1, NNod
        Id = (I - 1) * 2
        if (.not.IsFixDof(Id + 1).and.(MasS(Id + 1) .gt. 0.D0)) then
            !W(Id + 1) = W0(Id + 1) + (GrvFW(Id + 1) + ExtP(Id + 1) - InrFW(Id + 1) - DragF(Id + 1)) / MasW (Id + 1) * dtime
            temp=0.d0
            Temp = W0(Id + 1) + (GrvFW(Id + 1) + ExtP(Id + 1) - InrFW(Id + 1) - DragF(Id + 1) + DampW(Id + 1)) / MasW (Id + 1) * dtime
            W(Id + 1) = temp - sign(1.d0, temp) * dampf * abs(temp)
            aW = (W(Id + 1) - W0(Id + 1)) / dtime
            !V(Id + 1) = V0(Id + 1) + (GrvF(Id + 1) + ExtF(Id + 1) - InrF(Id + 1) - (aW * MasbW (Id + 1)) ) / MasS(Id + 1) * dtime
            temp=0.d0
            temp = V0(Id + 1) + (GrvF(Id + 1) + ExtF(Id + 1) - InrF(Id + 1) - (aW * MasbW (Id + 1)) + Damp(Id + 1) + DampW(Id + 1) ) / MasS(Id + 1) * dtime
            V(Id + 1) = temp - sign(1.d0, temp) * dampf * abs(temp)
        end if

        if (.not.IsFixDof(Id + 2).and.(MasS(Id + 2) .gt. 0.D0)) then
            !W(Id + 2) = W0(Id + 2) + (GrvFW(Id + 2) + ExtP(Id + 2) - InrFW(Id + 2) - DragF(Id + 2)) / MasW (Id + 2) * dtime
            temp = 0.d0
            Temp = W0(Id + 2) + (GrvFW(Id + 2) + ExtP(Id + 2) - InrFW(Id + 2) - DragF(Id + 2) + DampW(Id + 2)) / MasW (Id + 2) * dtime
            W(Id + 2) = temp - sign(1.d0, temp) * dampf * abs(temp)
            aW = (W(Id + 2) - W0(Id + 2)) / dtime
            !V(Id + 2) = V0(Id + 2) + (GrvF(Id + 2) + ExtF(Id + 2) - InrF(Id + 2) - (aW * MasbW (Id + 2)) ) / MasS(Id + 2) * dtime
            temp = 0.d0
            temp = V0(Id + 2) + (GrvF(Id + 2) + ExtF(Id + 2) - InrF(Id + 2) - (aW * MasbW (Id + 2)) + Damp(Id + 2) + DampW(Id + 2) ) / MasS(Id + 2) * dtime
            V(Id + 2) = temp - sign(1.d0, temp) * dampf * abs(temp)
        end if
    end do

    do IEl = 1, NEl
        INod(:) = ICon(:, IEl)
        do ig = 1, 4
            delV = 0.0
            delW = 0.0
            do I = 1, 4
                Id = (INod(I) - 1) * 2
                delV(1) = delV(1) + B(1, i, iel, ig) * V(Id + 1)
                delV(2) = delV(2) + B(2, i, iel, ig) * V(Id + 2)
                delV(3) = delV(3) + B(2, i, iel, ig) * V(Id + 1)
                delV(4) = delV(4) + B(1, i, iel, ig) * V(Id + 2)
                delW(1) = delW(1) + B(1, i, iel, ig) * W(Id + 1)
                delW(2) = delW(2) + B(2, i, iel, ig) * W(Id + 2)
                delW(3) = delW(3) + B(2, i, iel, ig) * W(Id + 1)
                delW(4) = delW(4) + B(1, i, iel, ig) * W(Id + 2)
            end do
            ! Strain increment
            dEps(1:2) = delV(1:2) * dtime
            dEps(3) = (delV(3) + delV(4)) * dtime

            ! Total strain
            EpsG(1:2, ig, IEl) = EpsG(1:2, ig, IEl) +  dEps(1:2)
            EpsG(3, ig, IEl) = EpsG(3, ig, IEl) +  dEps(3)

            ! Incremental strain
            dEpsG(1:2, ig, IEl) = dEps(1:2)
            dEpsG(3, ig, IEl) = dEps(3)

            ! Pressure increment
            dPore = (( (1-poro) * delV (1) + poro * delW(1) ) + ((1-poro) * delV (2) + poro * delW(2) ) ) * ( ( dtime * BulkW ) / poro )

            ! Total pressure
            PoreP(1:2, 1:4, IEl) = PoreP(1:2, 1:4, IEl) +  dPore
            PoreP(3, 1:4, IEl) = 0.d0

        enddo !gauss

        !do ig = 1, 4 ! gauss
        !    do I = 1, 4
        !        Id = (INod(I) - 1) * 2
        !        if(MasS(Id + 1).gt.0.d0) dEpsG(1:3, ig, IEl) = dEpsG(1:3, ig, IEl)*MasS(Id+1)
        !    end do
        !enddo !
        !
        !do ig = 1, 4 ! gauss
        !    do I = 1, 4
        !        Id = (INod(I) - 1) * 2
        !        if(MasS(Id + 1).gt.0.d0) dEpsG(1:3, ig, IEl) = dEpsG(1:3, ig, IEl)/MasS(Id+1)
        !    end do
        !enddo !gauss


        do ig = 1, 4
            
            call Elastic(MatProp(1), MatProp(2), dEpsg(:,ig,IEL), SigE(:,ig,IEl))
            !dSig = 0.d0
            !call Elastic(MatProp(1), MatProp(2), deps, dSig)
            !SigE(:,ig,IEl) = SigE(:,ig,IEl) + dSig(1:3)
            ! Plastic corrector, Mohr Coulomb. Using the Plaxis implementation here
            if(itime*dt.gt.5000.d0) then
                call MOHRC(SIN(4.d0*atan(1.d0)*10.d0/180.D0),0.d0,0.d0,MatProp(1),MatProp(2),SigE(:,ig,IEl),NPLAST,fbarnew(1,ig,iel),PlastInd(1,ig,Iel)) ! Mohr model
            end if
        enddo !gauss


        !do ig = 1, 4 ! gauss
        !    do I = 1, 4
        !        Id = (INod(I) - 1) * 2
        !        if(MasS(Id + 1).gt.0.d0) SigE(1:3, ig, IEl) = SigE(1:3, ig, IEl)*MasS(Id+1)
        !    end do
        !enddo !
        !
        !do ig = 1, 4 ! gauss
        !    do I = 1, 4
        !        Id = (INod(I) - 1) * 2
        !        if(MasS(Id + 1).gt.0.d0) SigE(1:3, ig, IEl) = SigE(1:3, ig, IEl)/MasS(Id+1)
        !    end do
        !enddo !gauss

        !SigE(:,ig,IEl) = SigE(:,ig,IEl) + dSig(1:3)
        do ig = 1, 4

            Sigg(:,ig,IEL) = SigE(:,ig,IEL) + PoreP(:,ig,IEL)
        enddo !gauss
        end do


        dis = dis + v * dt
        v0 = v
        w0 = w

    end subroutine Update_small

    !*********************************************************
    !Elastic Hookes model
    !**********************************************************
    subroutine Elastic(E, nu, eps, Sig)

    implicit none
    !
    double precision, intent(in) :: E, nu
    double precision, intent(inout) :: Sig(8)
    !double precision, intent(out) :: Sig(4)
    double precision, intent(in) :: eps(3)
    ! local variables
    double precision :: G_mod, K_mod, Eps_tr

    G_mod = E / (2.d0 * (1 + nu))
    K_mod = E / (3.d0 * (1 - 2 * nu))

    Eps_tr = eps(1) + eps(2)
    Sig(1) = ((K_mod * Eps_tr) + 2 * G_mod * (eps(1) - (Eps_tr/3.d0))) + Sig(1)
    Sig(2) = ((K_mod * Eps_tr) + 2 * G_mod * (eps(2) - (Eps_tr/3.d0))) + Sig(2)
    !Sig(3) = (2 * G_mod * eps(3))                                      + Sig(3)
    !Sig(4) = ((K_mod * Eps_tr) + 2 * G_mod * (0.d0 - (Eps_tr/3.d0)))   + Sig(4)

    Sig(3) = sig(3) +  ((K_mod * Eps_tr) + 2 * G_mod * (0.d0 - (Eps_tr/3.d0)))
    Sig(4) = sig(4) +  (2 * G_mod * eps(3))

    endsubroutine elastic



    subroutine VonMises(E, enu, Yield, Eps3, EpsP, s, EpsE)
    !*********************************************************************
    ! Von Mises Elasto-Plastic Model.
    !************************************************************************ ******

    implicit double precision (a - h, o - z)

    double precision, intent(in) :: E, enu, Yield, Eps3(3)
    double precision, intent(inout) :: S(7)
    double precision, intent(inout) :: EpsP(6), EpsE(3)
    integer :: i
    double precision :: Eps(6), devt(6), dir(6), dev(6), eps_v(3), EpsPT(3)

    B_K = E / (3.d0 * (1 - 2.d0 * enu)) ! Bulk-modulus K
    G_mod = E /(2.d0 * (1.d0 + enu)) ! 2nd lame - mu (G) shear mod

    gamma=0.d0
    eps = 0.d0
    eps(1:2) = Eps3(1:2)
    eps(4) = Eps3(3)

    ! write (LOgUnit, * ) Eps3
    treps = eps(1) + eps(2) + eps(3)

    ! dev of strain
    do i=1,3
        dev(i) = eps(i) - treps/3.d0
        dev(i+3) = eps(i+3) / 2.d0
    enddo

    !deviatoric trial force
    do i=1,6
        devt(i) = 2.d0 * G_mod * (dev(i) - EpsP(i))
    enddo

    !norm
    devtnorm = dsqrt(devt(1)**2.d0+ devt(2)**2.d0 + devt(3)**2.d0 + 2.d0* devt(4)**2.d0  &
        + 2.d0* devt(5)**2.d0 + 2.d0* devt(6)**2.d0)
    !write (LOGUnit, * ) devt(1), devt(2), devt(3)
    if (devtnorm.eq.0.d0)  devtnorm =0.0000001d0

    !direction of plastic flow
    do i = 1, 6
        dir(i) = devt(i)/(devtnorm)
    enddo

    !determine yield criterion
    Yn = (dsqrt(2.d0/3.d0) * (Yield ))
    phi = devtnorm - Yn

    !to compute stresses
    if ((phi .lt. 0.d0) ) then !elastic
        s(1) = B_K * treps + devt(1) !  stresses
        s(2) = B_K * treps + devt(2)
        s(3) = B_K * treps + devt(3)
        s(4) =  devt(4) !   plastic strains
        s(5) = EpsP(2)
        s(6) = EpsP(4)
        s(7) = gamma
        EpsE(1) = eps(1)
        EpsE(2) = eps(2)
        EpsE(3) = eps(3)

    else !plastic
        flag = 1.d0
        Yn = (dsqrt(2.d0/3.d0) * (Yield ))

        gamma = phi / (2.d0  * G_mod)

        !update plastic strain
        do i = 1, 6
            EpsP(i) = EpsP(i)+ (dir(i)*(gamma))
        enddo


        EpsPv = (Eps(1) - EpsE(1) + Eps(2) - EpsE(2) + Eps(3) - EpsE(3)) / 3.d0

        do i=1,3
            EpsPT(i) = EpsP(i) + EpsPv
        enddo

        treps1 =Eps(1) - EpsPT(1) + Eps(2) - EpsPT(2) + Eps(3) - EpsPT(3)

        s(1) = B_K * treps1 + devt(1) - ( 2.d0 * G_mod * flag * gamma * dir(1)) !  stresses
        s(2) = B_K * treps1 + devt(2) - ( 2.d0 * G_mod *  flag * gamma * dir(2))
        s(3) =  B_K * treps1 + devt(3) - ( 2.d0 * G_mod *  flag * gamma * dir(3))
        s(4) = devt(4) - 2.d0 * G_mod * (flag * gamma * dir(4)) !   plastic strains
        s(5) = EpsP(2)
        s(6) = EpsP(4)
        sigeq = (0.5d0 * ((s(1)-s(2))**2 + (s(2)-s(3))**2  +  (s(3)-s(1))**2 + 6*s(4)**2))
        s(7) = dsqrt (sigeq)

    endif

    end subroutine VonMises



    subroutine Eigen(a, x, abserr, n)
    !===========================================================
    ! Evaluate eigenvalues and eigenvectors
    ! of a real symmetric matrix a(n,n): a*x = lambda*x
    ! method: Jacoby method for symmetric matrices
    ! Alex G. (December 2009)
    !-----------------------------------------------------------
    ! input ...
    ! a(n,n) - array of coefficients for matrix A
    ! n      - number of equations
    ! abserr - abs tolerance [sum of (off-diagonal elements)^2]
    ! output ...
    ! a(i,i) - eigenvalues
    ! x(i,j) - eigenvectors
    ! comments ...
    !===========================================================
    implicit none
    integer i, j, k, n
    double precision a(n, n), x(n, n)
    double precision abserr, b2, bar
    double precision beta, coeff, c, s, cs, sc

    ! initialize x(i,j)=0, x(i,i)=1
    ! *** the array operation x=0.0 is specific for Fortran 90/95
    x = 0.0
    do i = 1, n
        x(i, i) = 1.0
    end do

    ! find the sum of all off-diagonal elements (squared)
    b2 = 0.0
    do i = 1, n
        do j = 1, n
            if (i .ne. j) b2 = b2 + a(i, j)**2
        end do
    end do

    if (b2 <= abserr) return

    ! average for off-diagonal elements /2
    bar = 0.5 * b2/float(n * n)

    do while (b2 .gt. abserr)
        do i = 1, n - 1
            do j = i + 1, n
                if (a(j, i)**2 <= bar) cycle ! do not touch small elements
                b2 = b2 - 2.0 * a(j, i)**2
                bar = 0.5 * b2/float(n * n)
                ! calculate coefficient c and s for Givens matrix
                beta = (a(j, j) - a(i, i))/(2.0 * a(j, i))
                coeff = 0.5 * beta/sqrt(1.0 + beta**2)
                s = sqrt(max(0.5 + coeff, 0.0))
                c = sqrt(max(0.5 - coeff, 0.0))
                ! recalculate rows i and j
                do k = 1, n
                    cs = c * a(i, k) + s * a(j, k)
                    sc = -s * a(i, k) + c * a(j, k)
                    a(i, k) = cs
                    a(j, k) = sc
                end do
                ! new matrix a_{k+1} from a_{k}, and eigenvectors
                do k = 1, n
                    cs = c * a(k, i) + s * a(k, j)
                    sc = -s * a(k, i) + c * a(k, j)
                    a(k, i) = cs
                    a(k, j) = sc
                    cs = c * x(k, i) + s * x(k, j)
                    sc = -s * x(k, i) + c * x(k, j)
                    x(k, i) = cs
                    x(k, j) = sc
                end do
            end do
        end do
    end do
    return
    end subroutine Eigen

    subroutine MkOpFiles()
    !**********************************************************************
    !    Function: Output files
    !**********************************************************************
    implicit none

    call MkFile(MSHUnit, 'Model.post.msh')
    call MkFile(RESUnit, 'Model.post.res')

    end subroutine MkOpFiles


    subroutine MkFile(Unit, flNam)
    !**********************************************************************
    !
    !    Function: Make a file
    !
    !**********************************************************************

    implicit none

    integer Unit
    character flNam * (*)

    if (FlExist(flNam)) then
        open(Unit, file = flNam)
        close(Unit, Status = 'Delete')
    endif
    open(Unit, file = flNam)

    end subroutine MkFile

    logical function FlExist(flNam)
    !**********************************************************************
    !
    !    Function: To check the existence of a file
    !
    !**********************************************************************

    implicit none

    logical lfil
    character flNam * (*)

    lfil = .false.
    inquire(file = flNam, exist = lfil)
    if (lfil) then
        FlExist = .true.
    else
        FlExist = .false.
    endif

    end function FlExist


    subroutine WtMSH(Unit)
    !**********************************************************************
    !    Function: Writing GiD *.msh file
    !!**********************************************************************

    implicit none
    integer Unit
    ! local variables
    integer :: IPart, INod, IEl, J, K
    double precision, dimension(2) :: Pos, r1, r2, VPos

    write(Unit, *) 'MESH dimension 2 ElemType Quadrilateral Nnode 4'
    write(Unit, *) 'Coordinates'
    do INod = 1, NNod
        write(Unit, 111) INod, NodCo(1, INod), NodCo(2, INod)
    end do
    write(Unit, *) 'End Coordinates'
    write(Unit, *) 'Elements'
    do IEl = 1, NEl
        write(Unit, "(5I7)") IEl, ICon(1, IEl), ICon(2, IEl), ICon(3, IEl), ICon(4, IEl)
    end do
    write(Unit, *) 'End Elements'
    close(Unit)
111 format(I7, 2E16.6E3)

    end subroutine WtMSH

    subroutine WtRES(IStep)
    !**********************************************************************
    !
    !    Function: Output to GiD
    !
    !**********************************************************************

    implicit none

    integer IStep
    ! local variables
    integer :: IEl, J, K, INod, Id, Unit, ig

    Unit = ResUnit
    if (IStep .eq. 1) then
        write(Unit, *) 'GiD Post Results File 1.0'
        write(Unit, *) 'GaussPoints "Material_Point" Elemtype Quadrilateral'
        write(Unit, *) 'Number of Gauss Points: 4'
        write(Unit, *) 'Natural Coordinates: Internal'
        write(Unit, *) 'end gausspoints'
        write(Unit, *) 'Result "Boundary" "MPM"', IStep, 'Vector OnNodes'
        write(Unit, *) 'ComponentNames "X-fix", "Y-fix"'
        write(Unit, *) 'values'
        do INod = 1, NNod
            Id = (INod - 1) * 2
            J = 0; K = 0
            if (IsFixDof(Id + 1)) J = 1
            if (IsFixDof(Id + 2)) K = 1
            write(Unit, "(3I7)") INod, J, K
        end do
        write(Unit, *) 'end values'
    end if
    write(Unit, *) 'Result "displacement" "MPM"', IStep, 'Vector OnNodes'
    write(Unit, *) 'ComponentNames "comp. x", "comp. y"'
    write(Unit, *) 'values'
    do INod = 1, NNod
        Id = (INod - 1) * 2
        write(Unit, "(I7, 2E16.6E3)") INod, Dis(Id + 1), Dis(Id + 2)
    end do
    write(Unit, *) 'end values'

    write(Unit, *) 'Result "Stress" "MPM"', IStep, 'Vector OnGaussPoints "Material_Point"'
    write(Unit, *) 'ComponentNames "sigma xx", "sigma yy", "sigma xy"'
    write(Unit, *) 'values'
    do IEl = 1, NEl
        write(Unit, "(I7, 4E16.6E4)") IEl, Sigg(1, 1, IEl), Sigg(2, 1, IEl), Sigg(3, 1, IEl)
        write(Unit, "(4E16.6E4)") Sigg(1, 2, IEl), Sigg(2, 2, IEl), Sigg(3, 2, IEl)
        write(Unit, "(4E16.6E4)") Sigg(1, 3, IEl), Sigg(2, 3, IEl), Sigg(3, 3, IEl)
        write(Unit, "(4E16.6E4)") Sigg(1, 4, IEl), Sigg(2, 4, IEl), Sigg(3, 4, IEl)
    end do
    write(Unit, *) 'end values'

    write(Unit, *) 'Result "EffStr" "MPM"', IStep, 'Vector OnGaussPoints "Material_Point"'
    write(Unit, *) 'ComponentNames "sigmae xx", "sigmae yy", "sigmae xy"'
    write(Unit, *) 'values'
    do IEl = 1, NEl
        write(Unit, "(I7, 4E16.6E4)") IEl, SigE(1, 1, IEl), SigE(2, 1, IEl), SigE(3, 1, IEl)
        write(Unit, "(4E16.6E4)") SigE(1, 2, IEl), SigE(2, 2, IEl), SigE(3, 2, IEl)
        write(Unit, "(4E16.6E4)") SigE(1, 3, IEl), SigE(2, 3, IEl), SigE(3, 3, IEl)
        write(Unit, "(4E16.6E4)") SigE(1, 4, IEl), SigE(2, 4, IEl), SigE(3, 4, IEl)
    end do
    write(Unit, *) 'end values'

    write(Unit, *) 'Result "Strain" "MPM"', IStep, 'Vector OnGaussPoints "Material_Point"'
    write(Unit, *) 'ComponentNames "eps xx", "eps yy", "eps xy"'
    write(Unit, *) 'values'
    do IEl = 1, NEl
        write(Unit, "(I7, 3E16.6E3)") IEl, EpsG(1, 1, IEl), EpsG(2, 1, IEl), EpsG(3, 1, IEl)
        write(Unit, "( 3E16.6E3)") Epsg(1, 2, IEl), EpsG(2, 2, IEl), EpsG(3, 2, IEl)
        write(Unit, "( 3E16.6E3)") Epsg(1, 3, IEl), Epsg(2, 3, IEl), Epsg(3, 3, IEl)
        write(Unit, "( 3E16.6E3)") Epsg(1, 4, IEl), Epsg(2, 4, IEl), Epsg(3, 4, IEl)
    end do
    write(Unit, *) 'end values'

    write(Unit, *) 'Result "Pore" "MPM"', IStep, 'Vector OnGaussPoints "Material_Point"'
    write(Unit, *) 'ComponentNames "Pore xx", "Pore yy", "Pore xy"'
    write(Unit, *) 'values'
    do IEl = 1, NEl
        write(Unit, "(I7, 4E16.6E4)") IEl, PoreP(1, 1, IEl), PoreP(2, 1, IEl), PoreP(3, 1, IEl)
        write(Unit, "(4E16.6E4)") PoreP(1, 2, IEl), PoreP(2, 2, IEl), PoreP(3, 2, IEl)
        write(Unit, "(4E16.6E4)") PoreP(1, 3, IEl), PoreP(2, 3, IEl), PoreP(3, 3, IEl)
        write(Unit, "(4E16.6E4)") PoreP(1, 4, IEl), PoreP(2, 4, IEl), PoreP(3, 4, IEl)
    end do
    write(Unit, *) 'end values'

    write(Unit, *) 'Result "PlastInd" "MPM"', IStep, 'Vector OnGaussPoints "Material_Point"'
    write(Unit, *) 'ComponentNames "PlastInd xx", "PlastInd yy", "PlastInd xy"'
    write(Unit, *) 'values'
    do IEl = 1, NEl
        write(Unit, "(I7, 4E16.6E4)") IEl, PlastInd(1, 1, IEl), PlastInd(2, 1, IEl), PlastInd(3, 1, IEl)
        write(Unit, "(4E16.6E4)") PlastInd(1, 2, IEl), PlastInd(2, 2, IEl), PlastInd(3, 2, IEl)
        write(Unit, "(4E16.6E4)") PlastInd(1, 3, IEl), PlastInd(2, 3, IEl), PlastInd(3, 3, IEl)
        write(Unit, "(4E16.6E4)") PlastInd(1, 4, IEl), PlastInd(2, 4, IEl), PlastInd(3, 4, IEl)
    end do
    write(Unit, *) 'end values'

    end subroutine WtRES

    end module modulefem

    !************program *****************
    program MPM2D

    use modulefem
    implicit none

    call MkOpFiles()
    call readdata()

    call solve()

    end program
    !*************program*****************
