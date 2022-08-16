module array
    implicit none
    double precision, allocatable,dimension(:)::s_delta,delta,mi,th,ph,rwork,evl_s,n_i,s_x,s_y,s_z,n_i_up,n_i_dn
    integer,dimension(:),allocatable::x(:),y(:),tag(:)
    complex*16,allocatable:: h(:,:),bb(:,:),work(:),mat(:,:),h_(:,:)
    double precision::nsum_t(1600),nn1(1600),ee1(1600),e_2(1600)
    doubleprecision,allocatable::phi(:),theta(:)
    complex*16,allocatable::mat_1(:,:),mat_2(:,:),mat_3(:,:),vec(:) 
    complex*16,allocatable::mat_4(:,:),mat_5(:,:),mat_6(:,:),alpha(:)

	!complex*16,allocatable,dimension(:,:):: h3,h1,h4,h2
    complex*16,allocatable::mat_1_up(:,:),mat_2_up(:,:),mat_3_up(:,:)
    complex*16,allocatable::mat_4_up(:,:),mat_5_up(:,:),mat_6_up(:,:)
    complex*16,allocatable::mat_4_down(:,:),mat_5_down(:,:),mat_6_down(:,:)
    complex*16,allocatable::mat_1_down(:,:),mat_2_down(:,:),mat_3_down(:,:),mat_temp(:,:),mat_temp_down(:,:)
 
end module array

module s_vd
    complex*16,allocatable,dimension(:,:)::A1,U_1,VT
    !complex*16,allocatable,dimension(:,:)::
    complex*16, allocatable::sigma1(:)
endmodule s_vd

module global  
        implicit none  
        double precision::t,ua,ub,uc,mu,m_sum,e_sum,u,esum_t,m1,ph1,th1,nsum,e,s_f,s_f1,s_f2
        double precision::lambda,sx_i_sx_j,sy_i_sy_j,sz_i_sz_j,ch_no,lambda_r
        integer::nos,unit_cells
        integer::i_1,i_2,i_3,i,j,k_x,k_y,b1,d,ie,sigma,a1,lx1,ly1,jx,jy,lx,ly,nx,ny,dkx,dky
endmodule global

program liebmain
    use array
    use global
    use mtmod
    implicit none
    integer::seed,count1,flag1
    double precision::temp,get_mu_s,de,ep,ec,s,bott_index
    double precision::inpt,filling,fill_inpt

    open(8,file='input.dat',status='unknown')
    do i=1,8
    read(8,*)inpt
    if (i.eq.1) seed=int(inpt)
    if (i.eq.2) d=int(inpt)
    if (i.eq.3) t=dble(inpt)
    if (i.eq.4) ua=dble(inpt)
    if (i.eq.5) ub=dble(inpt)
    if (i.eq.6) uc=dble(inpt)
    if (i.eq.7) fill_inpt=dble(inpt)
    if (i.eq.8)  lambda=dble(inpt)
    end do
    close(8)
    lambda_r=0.30d0
    call sgrnd(seed)

    unit_cells=(d/2)**2
    nos=3*unit_cells
    filling=fill_inpt*2*nos    

    allocate(x(nos),y(nos),h((2*nos),(2*nos)),tag(nos))
    allocate(mat_temp((d/2)**2,(2*nos)),mat_temp_down((d/2)**2,(2*nos)))
    allocate(evl_s(2*nos))
    allocate(n_i((nos)),delta(nos),s_delta(nos))
    allocate(s_x((nos)))
    allocate(s_y((nos)))
    allocate(s_z((nos)))
    allocate(n_i_up((nos)))
    allocate(n_i_dn((nos)))
    allocate(phi(nos),theta(nos))
    allocate(alpha(2*nos))

    print*,'unit cells_____________=',unit_cells
    print*,'sites__________________=',nos
    print*,'Matrix dim_____________=',2*nos,'X',2*nos
    print*,'hopping parameter______=',t
    print*,'Ua_____________________=',ua
    print*,'ub_____________________=',ub
    print*,'uc_____________________=',uc
    print*,'filling________________=',fill_inpt,filling
    print*,'lambda_________________=',lambda

    temp=1e-04
    call lattice_labeling
    call stagg_v
    call matgen
    !call matgen_so
    call matgen_hf_u
    call matgen_rsoc

    call check_Hermiticity
    call diagonalization(h,2)
    do i=1,2*nos
        write(726,*) i, evl_s(i)
    enddo
    mu=get_mu_s(filling,temp) 
    call dos
    call bloch_states6(h)
    do i=1,6
    call chern_input
    call chern_number
    enddo

    do i_1=1,6
        call bottindex(i_1,i_1+1,bott_index)
        print*,i_1,'band','________BI=',bott_index
    enddo
    stop
    ! call bloch_states6(h)
    ! mat_temp=mat_1
    ! call chern_number
    ! stop
    ! call bloch_states(h)
    ! mat_temp=mat_1_down
    ! call chern_number
    ! stop

    call lattice_labeling
    call initial_parameter
    do i=1,nos
        write(86,*)i,n_i(i),s_x(i),s_y(i),s_z(i),delta(i)
    enddo       
    call matgen
    call matgen_so
    call matgen_hf_u
    call diagonalization(h,2)
    mu=get_mu_s(filling,temp)  

    call output_parameter(temp) 
    call energy_cal(temp)
    Ep=esum_t
    count1=0
    flag1=0
        
    do while(flag1.eq.0)
        
        count1=count1+1
        call matgen
        call matgen_so
        call matgen_hf_u
        call diagonalization(h,2)
        mu=get_mu_s(filling,temp)
        call energy_cal(temp)
        Ec=esum_t
        de=Ec-Ep
        write(94,*)count1,ec,ep
        flush(94)

        if(abs(de).lt.1e-10)then
            flag1=1
            call output_parameter(temp)
            do i=1,nos
                write(500,*)i,n_i(i),s_x(i),s_y(i),s_z(i),n_i_up(i),n_i_dn(i)
                write(501,*)i,sqrt(s_x(i)*s_x(i)+s_y(i)*s_y(i)+s_z(i)*s_z(i))
            enddo           
            call matgen
            call matgen_so
            call matgen_hf_u
            call diagonalization(h,2)
            ! do i=1,2*nos
            !     write(726,*) i, evl_s(i)
            ! enddo
            mu=get_mu_s(filling,temp)  
            
            do i_1=1,6
                call bottindex(i_1,i_1+1,bott_index)
                print*,i_1,'band','________BI=',bott_index
            enddo

            ! call bloch_states6(h)
            ! do i=1,6
            !     call chern_input
            !     call chern_number
            ! enddo

            ! call bottindex(1,2,bott_index)
            ! print*,bott_index

            
            ! call bloch_states6(h)
            ! mat_temp=mat_1_down
            ! call chern_number   
                    
            call dos
            ! call proj_dos(temp)
            ! call occupation(temp)
            ! call s_xy_corr_a(temp)
            ! call sz_corr_a(temp)
            
            ! s=sqrt(sx_i_sx_j+sy_i_sy_j+sz_i_sz_j)
            ! write(181,*) s
        else
            flag1=0
            Ep=esum_t
            call output_parameter(temp)
        endif            
    enddo 

    write(61,*) count1
end program liebmain

subroutine chern_input
    use global
    use array
    implicit none
  
  
    if(i.eq.1) mat_temp=mat_1    
    if(i.eq.2) mat_temp=mat_2  
    if(i.eq.3) mat_temp=mat_3
    if(i.eq.4) mat_temp=mat_4
    if(i.eq.5) mat_temp=mat_5
    if(i.eq.6) mat_temp=mat_6
    ! if(i.eq.7) mat_temp=mat_4_up     
    ! if(i.eq.8) mat_temp=mat_4_down
    ! if(i.eq.9) mat_temp=mat_5_up    
    ! if(i.eq.10) mat_temp=mat_5_down   
    ! if(i.eq.11) mat_temp=mat_6_up    
    ! if(i.eq.12) mat_temp=mat_6_down
  
  
endsubroutine chern_input

subroutine lattice_labeling
    use global
    use array
    implicit none
    integer::ix,iy,i1
    double precision::Pi
    pi=4.0*atan(1.0)
    i1=-2
    do iy=1,d
        do ix=1,d
            if((mod(iy,2).ne.0).and.(mod(ix,2).ne.0))then   

                i1=i1+3
                x(i1)=ix                     !A
                y(i1)=iy
                tag(i1)=1

                x(i1+1)=ix+1                 !B
                y(i1+1)=iy
                tag(i1+1)=2

                x(i1+2)=ix                   !C
                y(i1+2)=iy+1
                tag(i1+2)=3

            endif
        enddo
    enddo
    i1=0
    do i=1,nos           
        if(tag(i).eq.1)then
            write(12,*) i,x(i),y(i),tag(i)
            flush(12)
        endif
        if(tag(i).eq.2)then 
 
            write(13,*) i,x(i),y(i),tag(i)
            flush(13)

        endif      
        if(tag(i).eq.3)then
   
            write(14,*) i,x(i),y(i),tag(i)
            flush(14) 
        endif     
    enddo

endsubroutine lattice_labeling


subroutine matgen
    use array
    use global
    implicit none
    integer::l,k,xi,yi,xd,yd,a,b

    h=complex(0.0d0,0.0d0)
    do l=1,nos
        if(tag(l).eq.1) then
            xi=1
            yi=1
            xd=1
            yd=1
            if(x(l).eq.d) xi=-d+1
            if(y(l).eq.d) yi=-d+1

            if(x(l).eq.1) xd=-d+1
            if(y(l).eq.1) yd=-d+1


            ! call rannum(delta_t)
            ! t=-1+3*delta_t

            do k=1,nos
                if((tag(k).eq.2).or.((tag(k).eq.3)))then        
                    if((x(k).eq.(x(l)+xi)).and.(y(k).eq.y(l)))then  
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=t;h(b,a)=t
                        h(a+1,b+1)=t;h(b+1,a+1)=t
                    endif

                    if((x(k).eq.(x(l)-xd)).and.(y(k).eq.y(l)))then  
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=t;h(b,a)=t
                        h(a+1,b+1)=t;h(b+1,a+1)=t
                    endif
        
                    if((x(k).eq.x(l)).and.(y(k).eq.(y(l)+yi)))then  
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=t;h(b,a)=t
                        h(a+1,b+1)=t;h(b+1,a+1)=t
                    endif
        
                    if((x(k).eq.x(l)).and.(y(k).eq.(y(l)-yd)))then  
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=t;h(b,a)=t
                        h(a+1,b+1)=t;h(b+1,a+1)=t
                    endif
  
                endif
            enddo
        endif
    enddo

 !____________zeeman field_______________________________

    ! do i=1,2*nos
    !     if(mod(i,2).ne.0)   h(i,i)=-0.5*0.8
    !     if(mod(i,2).eq.0)   h(i,i)=0.5*0.8
    ! enddo
 !______________________________________________________



endsubroutine matgen

subroutine check_Hermiticity
    use array
    use global
    implicit none
    integer::l,k
    do l=1,2*nos
        do k=1,2*nos
            !if(h(l,k).ne.0) write(34,*) l,k,h(l,k)
            if(h(l,k).ne.conjg(h(k,l))) print*,l,k,'not hermitian'
        enddo
    enddo 
endsubroutine check_Hermiticity

subroutine matgen_hf_u
    use array
    use global
    use mtmod
    implicit none
    integer::l,a,b

    do l=1,nos
        if(tag(l).eq.1) u=ua
        if(tag(l).eq.2) u=ub
        if(tag(l).eq.3) u=uc
        
        a=2*l-1
        b=2*l-1
        !print*,tag(l),delta(l),s_delta(l)
        h(a,b)=(u/2.0d0)*((-(s_z(l))*2.0d0)+n_i(l))+delta(l)+s_delta(l)
        h(a,b+1)=(u/2.0d0)*complex(((-s_x(l))*2.0d0),(s_y(l))*2.0d0)
        h(a+1,b)=conjg(h(a,b+1))
        h(a+1,b+1)=(u/2.0d0)*(((s_z(l))*2.0d0)+n_i(l))+delta(l)-s_delta(l)
        
    enddo

endsubroutine matgen_hf_u



subroutine matgen_so
    use global
    use array
    implicit none
    integer::l,a,b,k,xi,yi,xd,yd

    do l=1,nos
        if(tag(l).eq.2)then
            xi=1
            yi=1
            xd=1
            yd=1
            if(x(l).eq.d) xi=-d+1
            if(y(l).eq.d) yi=-d+1

            if(x(l).eq.1) xd=-d+1
            if(y(l).eq.1) yd=-d+1
            do k=1,nos
                if((tag(k).eq.3))then

                    if((x(k).eq.(x(l)+xi)).and.(y(k).eq.(y(l)+yi)))then 
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=-lambda*cmplx(0,1);h(b,a)=conjg(h(a,b))
                        h(a+1,b+1)=lambda*cmplx(0,1);h(b+1,a+1)=conjg(h(a+1,b+1))
                    endif

                    if((x(k).eq.(x(l)-xd)).and.(y(k).eq.(y(l)+yi)))then 
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=lambda*cmplx(0,1);h(b,a)=conjg(h(a,b))
                        h(a+1,b+1)=-lambda*cmplx(0,1);h(b+1,a+1)=conjg(h(a+1,b+1))

                    endif 

                    if((x(k).eq.(x(l)-xd)).and.(y(k).eq.(y(l)-yd)))then 
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=-lambda*cmplx(0,1);h(b,a)=conjg(h(a,b))
                        h(a+1,b+1)=lambda*cmplx(0,1);h(b+1,a+1)=conjg(h(a+1,b+1))
                    endif 

                    if((x(k).eq.(x(l)+xi)).and.(y(k).eq.(y(l)-yd)))then 
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=lambda*cmplx(0,1);h(b,a)=conjg(h(a,b))
                        h(a+1,b+1)=-lambda*cmplx(0,1);h(b+1,a+1)=conjg(h(a+1,b+1))

                    endif   

                endif
            enddo
        endif
    enddo
endsubroutine matgen_so

subroutine matgen_rsoc
    use array
    use global
    implicit none
    integer::l,k,xi,yi,xd,yd,a,b

    do l=1,nos
        if(tag(l).eq.1) then
            xi=1
            yi=1
            xd=1
            yd=1
            if(x(l).eq.d) xi=-d+1
            if(y(l).eq.d) yi=-d+1

            if(x(l).eq.1) xd=-d+1
            if(y(l).eq.1) yd=-d+1

            do k=1,nos
                if((tag(k).eq.2).or.((tag(k).eq.3)))then        
                    if((x(k).eq.(x(l)+xi)).and.(y(k).eq.y(l)))then  
                        a=2*l-1
                        b=2*k-1
                        h(a,b+1)=cmplx(0.0d0,1.0d0)*cmplx(0.0d0,1.0d0)*lambda_r;h(b+1,a)=conjg(h(a,b+1))
                        h(a+1,b)=-cmplx(0.0d0,1.0d0)*cmplx(0.0d0,1.0d0)*lambda_r;h(b,a+1)=conjg(h(a+1,b))
                    endif

                    if((x(k).eq.(x(l)-xd)).and.(y(k).eq.y(l)))then  
                        a=2*l-1
                        b=2*k-1
                        h(a,b+1)=-cmplx(0.0d0,1.0d0)*cmplx(0.0d0,1.0d0)*lambda_r;h(b+1,a)=conjg(h(a,b+1))
                        h(a+1,b)=cmplx(0.0d0,1.0d0)*cmplx(0.0d0,1.0d0)*lambda_r;h(b,a+1)=conjg(h(a+1,b))
                    endif
        
                    if((x(k).eq.x(l)).and.(y(k).eq.(y(l)+yi)))then  
                        a=2*l-1
                        b=2*k-1
                        h(a,b+1)=cmplx(0.0d0,1.0d0)*lambda_r;h(b+1,a)=conjg(h(a,b+1))
                        h(a+1,b)=cmplx(0.0d0,1.0d0)*lambda_r;h(b,a+1)=conjg(h(a+1,b))
                    endif
        
                    if((x(k).eq.x(l)).and.(y(k).eq.(y(l)-yd)))then  
                        a=2*l-1
                        b=2*k-1
                        h(a,b+1)=-cmplx(0.0d0,1.0d0)*lambda_r;h(b+1,a)=conjg(h(a,b+1))
                        h(a+1,b)=-cmplx(0.0d0,1.0d0)*lambda_r;h(b,a+1)=conjg(h(a+1,b))
                    endif
  
                endif
            enddo
        endif
    enddo

endsubroutine matgen_rsoc

subroutine diagonalization(h_temp,flag_diag)
    use array
    use global
    implicit none
    integer::lda,lwmax,info,lwork,flag_diag
    complex*16::h_temp(2*nos,2*nos)
    !EXTERNAL::SELECT
    allocate(rwork(3*(2*nos)-2),work(2*(2*nos)-1))
    lda=(2*nos)
    lwmax=(2*nos)
    lwork=(2*(2*nos)-1)
    evl_s=0.0d0
    if(flag_diag==1)then
       call zheev('n','u',2*nos,h_temp,lda,evl_s,work,lwork,rwork,info)
    endif

    if(flag_diag==2)then
       call zheev('v','u',2*nos,h_temp,lda,evl_s,work,lwork,rwork,info)
    endif

    if(info.ne.0)then
        print*,'algorithm failed'  
    endif 
   
    ! do i=1,2*nos
    !     write(726,*) i, evl_s(i)
    ! enddo
    
    deallocate(rwork,work)  
endsubroutine diagonalization

subroutine diagonalization1(h_temp)
    use array
    use global
    implicit none
    integer::lda,lwmax,info,lwork,ldvl,ldvr
    complex*16::h_temp(2*nos,2*nos),vs(2*nos),vl(2*nos,2*nos),vr(2*nos,2*nos)

    allocate(rwork(2*(2*nos)),work(2*(2*nos)))
    lda=(2*nos)
    lwmax=(2*nos)
    lwork=2*2*nos
    ldvl=2*nos
    ldvr=2*nos
    alpha=cmplx(0.0d0,0.0d0)
  
    lwork=2*(2*nos)
    call ZGEEV('n','n',2*nos,h_temp,LDA,alpha,VL,LDVL,VR,LDVR,WORK,LWORK,RWORK,INFO)
   
    if(info.ne.0)then
        print*,'algorithm failed'  
    endif 

    do i=1,2*nos
        write(727,*) i, alpha(i)
    enddo
    
    deallocate(rwork,work)  
endsubroutine diagonalization1

subroutine rannum(r)
    use mtmod
    implicit none
    double precision::r
    r=grnd()
    return
endsubroutine rannum

double precision function get_mu_s(fill,temp2)
    use array
    use global
    implicit none
    double precision::f,fL2,fR,mR,mL,m_d,temp2
    double precision::fill

    mR = maxval(evl_s)       !right-side chemical potential
    fr=0.0d0
    do i=1,(2*nos)
        fr=fr+(1.0d0/(exp((evl_s(i)-mR)/Temp2)+1.0d0))
    end do

    mL = minval(evl_s)       !left-side chemical potential
    fL2=0.0d0

    do i=1,(2*nos)
        fL2=fL2+(1.0d0/(exp((evl_s(i)-mL)/Temp2)+1.0d0))
    end do

    m_d = 0.5d0*(mL+mR)    !middle chemical potential
    f=0.0d0
    
    do i=1,(2*nos)
        f=f+(1.0d0/(exp((evl_s(i)-m_d)/Temp2)+1.0d0))
    end do     
    !print*,f,fill
    do while(abs(f-fill).ge.1e-8)
        m_d = 0.5d0*(mL+mR)
        f=0.0d0
        do i=1,(2*nos)
            f=f+(1.0d0/(exp((evl_s(i)-m_d)/Temp2)+1.0d0))
        end do
        if(f.gt.fill)then
            !if middle filling is above target, make it the new right bound.
            mR = m_d
            fR = f
        elseif(f.lt.fill)then
            !if middle filling is below target, make it the new left bound.
            mL = m_d
            fR = f
        endif
    enddo
        
    !Return the middle value
    get_mu_s = m_d
    return
end function get_mu_s

subroutine energy_cal(temp2)
    use array
    use global
    implicit none
    doubleprecision::temp2,esum
    esum=0.0d0

    do i=1,(2*nos)
        esum=esum+(evl_s(i)/(1.0d0+(exp((evl_s(i)-mu)/temp2))))
    enddo   
    esum_t=(esum)/(2*nos)       
endsubroutine energy_cal

subroutine output_parameter(temp2)
    use global
    use array
    implicit none
    doubleprecision::ni_up,ni_dn,s_x_sum,s_y_sum,fermi_fn,temp2

    s_x=0.0d0
    s_y=0.0d0
    n_i=0.0d0
    s_z=0.0d0
    n_i_up=0.0d0
    n_i_dn=0.0d0

    do i=1,nos

        ni_up=0.0d0
        ni_dn=0.0d0
        s_x_sum=0.0d0
        s_y_sum=0.0d0

        do j=1,(2*nos)
           
            fermi_fn=(1/(1.0d0+(exp((evl_s(j)-mu)/temp2))))

            ni_up=ni_up+(h(2*i-1,j)*conjg(h(2*i-1,j))*fermi_fn)
           
            ni_dn=ni_dn+(h(2*i,j)*conjg(h(2*i,j))*fermi_fn)
           
            s_x_sum=s_x_sum+((conjg(h(2*i-1,j))*h(2*i,j)+(h(2*i-1,j))*conjg(h(2*i,j)))*fermi_fn)
           
            s_y_sum=s_y_sum+(((conjg(h(2*i-1,j))*h(2*i,j)-(h(2*i-1,j))*conjg(h(2*i,j)))*fermi_fn)/(complex(0.0d0,1.0d0)))
        enddo
    
        ! s_x(i)=s_x_sum/2
        ! s_y(i)=s_y_sum/2
        ! s_z(i)=(ni_up-ni_dn)/2
        ! n_i(i)=ni_up+ni_dn
        ! n_i_up(i)=ni_up
        ! n_i_dn(i)=ni_dn

        s_x(i)=0
        s_y(i)=0
        n_i(i)=ni_up+ni_dn
        n_i_up(i)=ni_up
        n_i_dn(i)=ni_dn
        if(tag(i).eq.1) s_z(i)=0.5
        if(tag(i).eq.2) s_z(i)=0.5
        if(tag(i).eq.3) s_z(i)=0.5
    enddo      
endsubroutine output_parameter

subroutine initial_parameter
   use mtmod
   use global
   use array
   implicit none
   double precision::s_x_i,s_y_i,s_z_i,n_i_i
   double precision::v1

    do i=1,nos
        call rannum(n_i_i)
        call rannum(s_x_i)
        call rannum(s_y_i)
        call rannum(s_z_i)
        
        ! n_i(i)=n_i_i
        ! s_x(i)=s_x_i
        ! s_y(i)=s_y_i
        ! s_z(i)=s_z_i

        if(tag(i).eq.1)then
            if(mod(i,2).eq.0) delta(i)=0.0
            if(mod(i,2).ne.0) delta(i)=0.0
        endif

        if(tag(i).eq.2)then
            if(mod(i,2).eq.0) delta(i)=-1
            if(mod(i,2).ne.0) delta(i)=1
        endif

        if(tag(i).eq.3)then
            if(mod(i,2).eq.0) delta(i)=-1
            if(mod(i,2).ne.0) delta(i)=1
        endif
            
  

        n_i(i)=1
        s_x(i)=0
        s_y(i)=0
        if(tag(i).eq.1) s_z(i)=0.5
        if(tag(i).eq.2) s_z(i)=0.5
        if(tag(i).eq.3) s_z(i)=0.5
    enddo

endsubroutine initial_parameter

subroutine stagg_v
    use mtmod
    use global
    use array
    implicit none

    do i=1,nos

        if(tag(i).eq.1) delta(i)=0.2
        if(tag(i).eq.2) delta(i)=-0.2  
        if(tag(i).eq.3) delta(i)=-0.2

        if((tag(i).eq.1).and.(mod(i,2).eq.0)) s_delta(i)=0.1
        if((tag(i).eq.1).and.(mod(i,2).ne.0)) s_delta(i)=-0.1

        if((tag(i).eq.2).and.(mod(i,2).eq.0)) s_delta(i)=-0.1
        if((tag(i).eq.2).and.(mod(i,2).ne.0)) s_delta(i)=0.1

        if((tag(i).eq.3).and.(mod(i,2).eq.0)) s_delta(i)=-0.1
        if((tag(i).eq.3).and.(mod(i,2).ne.0)) s_delta(i)=0.1
    enddo
 
endsubroutine stagg_v

subroutine dos
    use array
    use global
    implicit none
    double precision::pi,eta_,wmin,wmax,dw,w
    integer::nwp

    pi=4.0*atan(1.0)
    wmin=-20
    wmax=20
    nwp=2000
    dw=abs(wmax-wmin)/nwp
    w=wmin
    eta_=0.10d0

    do ie=1,nwp!no of interval bw bandwitdh of e	
        w=w+dw
        nsum=0.0d0

        do j=1,2*nos
           nsum=nsum+((eta_/pi)/((w-(evl_s(j)-mu))**2+((eta_)**2)))
        enddo
   
        write(200,*) w,nsum/(2*nos)
        flush(200)
    enddo  
           
endsubroutine dos

subroutine s_xy_corr_b_c(temp2)
    use global
    use array
    implicit none
    integer::j1,l
    doubleprecision::c_i_dn_cdag_i_up,c_j_dn_cdag_j_up,c_i_dn_cdag_j_up,c_i_up_cdag_i_dn,cdag_i_dn_c_j_up 
    double precision::cdag_i_up_c_j_up,cdag_i_dn_c_j_dn,cdag_i_up_c_j_dn,c_i_up_cdag_j_dn,c_j_up_cdag_j_dn
    doubleprecision::c_i_up_cdag_j_up,fermi_fn,sx1,sx2,sx3,sx4,temp2

    sx_i_sx_j=0.0d0
    sy_i_sy_j=0.0d0
    do i=1,nos
    do j=1,nos
        if(i==j)then
        if((mod((x(i)+y(i)),2).ne.0).and.(mod((x(j)+y(j)),2).ne.0))then  ! b & c 
                    
            
        c_i_dn_cdag_i_up=0.0d0
        c_j_dn_cdag_j_up=0.0d0
        c_i_dn_cdag_j_up=0.0d0
        c_i_up_cdag_i_dn=0.0d0
        c_i_up_cdag_j_up=0.0d0
        c_j_up_cdag_j_dn=0.0d0
        c_i_up_cdag_j_dn=0.0d0
        cdag_i_up_c_j_dn=0.0d0
        cdag_i_dn_c_j_dn=0.0d0
        cdag_i_up_c_j_up=0.0d0
        cdag_i_dn_c_j_up=0.0d0


        do l=1,2*nos

            fermi_fn=(1/(1.0d0+(exp((evl_s(l)-mu)/temp2))))

            c_i_dn_cdag_i_up=c_i_dn_cdag_i_up+((h(2*i,l))*conjg(h(2*i-1,l)))*(1-fermi_fn)   !1

            c_j_dn_cdag_j_up=c_j_dn_cdag_j_up+((h(2*j,l))*conjg(h(2*j-1,l)))*(1-fermi_fn)  !2

            c_i_dn_cdag_j_up=c_i_dn_cdag_j_up+((h(2*i,l))*conjg(h(2*j-1,l)))*(1-fermi_fn)   !3

            c_i_up_cdag_i_dn=c_i_up_cdag_i_dn+((h(2*i-1,l))*conjg(h(2*i,l)))*fermi_fn    !4

            c_i_up_cdag_j_up=c_i_up_cdag_j_up+((h(2*i-1,l))*conjg(h(2*j-1,l)))*fermi_fn   !6

            c_j_up_cdag_j_dn=c_j_up_cdag_j_dn+((h(2*j-1,l))*conjg(h(2*j,l)))*fermi_fn     !9

            c_i_up_cdag_j_dn=c_i_up_cdag_j_dn+((h(2*i-1,l))*conjg(h(2*j,l)))*fermi_fn     !10

            cdag_i_up_c_j_dn=cdag_i_up_c_j_dn+((h(2*j,l))*conjg(h(2*i-1,l)))*fermi_fn      !5

            cdag_i_dn_c_j_dn=cdag_i_dn_c_j_dn+((h(2*j,l))*conjg(h(2*i,l)))*fermi_fn       !7

            cdag_i_up_c_j_up=cdag_i_up_c_j_up+((h(2*j-1,l))*conjg(h(2*i-1,l)))*fermi_fn   !8

            cdag_i_dn_c_j_up=cdag_i_dn_c_j_up+((h(2*j-1,l))*conjg(h(2*i,l)))*fermi_fn     !11

        enddo

        
        
        sx1=(c_i_dn_cdag_i_up*c_j_dn_cdag_j_up)+(c_i_dn_cdag_j_up*cdag_i_up_c_j_dn)

        sx2=(c_i_up_cdag_i_dn*c_j_dn_cdag_j_up)+(c_i_up_cdag_j_up*cdag_i_dn_c_j_dn)

        sx3=(c_i_dn_cdag_i_up*c_j_up_cdag_j_dn)+(c_i_dn_cdag_j_up*cdag_i_up_c_j_up)

        sx4=(c_i_up_cdag_i_dn*c_j_up_cdag_j_dn)+(c_i_up_cdag_j_dn*cdag_i_dn_c_j_up)

        sx_i_sx_j=sx_i_sx_j+(sx1+sx2+sx3+sx4)*((-1)**(abs(i-j)))

        sy_i_sy_j=sy_i_sy_j+(-sx1+sx2+sx3-sx4)*((-1)**(abs(i-j)))

           endif
           endif
    enddo
enddo

!print*,'Sx_Sx_b_c',sx_i_sx_j !/(nos/2)
!print*,'Sy_Sy_b_c',sy_i_sy_j !/(nos/3)

        
endsubroutine s_xy_corr_b_c	 


subroutine s_xy_corr_a(temp2)
    use global
    use array
    implicit none
    integer::j1,l
    doubleprecision::c_i_dn_cdag_i_up,c_j_dn_cdag_j_up,c_i_dn_cdag_j_up,c_i_up_cdag_i_dn,cdag_i_dn_c_j_up 
    double precision::cdag_i_up_c_j_up,cdag_i_dn_c_j_dn,cdag_i_up_c_j_dn,c_i_up_cdag_j_dn,c_j_up_cdag_j_dn
    doubleprecision::c_i_up_cdag_j_up,fermi_fn,sx1,sx2,sx3,sx4,temp2

    sx_i_sx_j=0.0d0
    sy_i_sy_j=0.0d0
    
    do i=1,nos
        if(tag(i).eq.1)then
            do j=1,nos
                if((tag(j).eq.1).and.(i.eq.j))then
                    c_i_dn_cdag_i_up=0.0d0
                    c_j_dn_cdag_j_up=0.0d0
                    c_i_dn_cdag_j_up=0.0d0
                    c_i_up_cdag_i_dn=0.0d0
                    c_i_up_cdag_j_up=0.0d0
                    c_j_up_cdag_j_dn=0.0d0
                    c_i_up_cdag_j_dn=0.0d0
                    cdag_i_up_c_j_dn=0.0d0
                    cdag_i_dn_c_j_dn=0.0d0
                    cdag_i_up_c_j_up=0.0d0
                    cdag_i_dn_c_j_up=0.0d0


                    do l=1,2*nos

                        fermi_fn=(1/(1.0d0+(exp((evl_s(l)-mu)/temp2))))

                        c_i_dn_cdag_i_up=c_i_dn_cdag_i_up+((h(2*i,l))*conjg(h(2*i-1,l)))*(1-fermi_fn)   !1

                        c_j_dn_cdag_j_up=c_j_dn_cdag_j_up+((h(2*j,l))*conjg(h(2*j-1,l)))*(1-fermi_fn)    !2

                        c_i_dn_cdag_j_up=c_i_dn_cdag_j_up+((h(2*i,l))*conjg(h(2*j-1,l)))*(1-fermi_fn)    !3

                        c_i_up_cdag_i_dn=c_i_up_cdag_i_dn+((h(2*i-1,l))*conjg(h(2*i,l)))*(1-fermi_fn)    !4

                        c_i_up_cdag_j_up=c_i_up_cdag_j_up+((h(2*i-1,l))*conjg(h(2*j-1,l)))*(1-fermi_fn)   !6

                        c_j_up_cdag_j_dn=c_j_up_cdag_j_dn+((h(2*j-1,l))*conjg(h(2*j,l)))*(1-fermi_fn)     !9

                        c_i_up_cdag_j_dn=c_i_up_cdag_j_dn+((h(2*i-1,l))*conjg(h(2*j,l)))*(1-fermi_fn)     !10

                        cdag_i_up_c_j_dn=cdag_i_up_c_j_dn+((h(2*j,l))*conjg(h(2*i-1,l)))*fermi_fn      !5

                        cdag_i_dn_c_j_dn=cdag_i_dn_c_j_dn+((h(2*j,l))*conjg(h(2*i,l)))*fermi_fn       !7

                        cdag_i_up_c_j_up=cdag_i_up_c_j_up+((h(2*j-1,l))*conjg(h(2*i-1,l)))*fermi_fn   !8

                        cdag_i_dn_c_j_up=cdag_i_dn_c_j_up+((h(2*j-1,l))*conjg(h(2*i,l)))*fermi_fn     !11

                    enddo

        
                    sx1=(c_i_dn_cdag_i_up*c_j_dn_cdag_j_up)+(c_i_dn_cdag_j_up*cdag_i_up_c_j_dn)

                    sx2=(c_i_up_cdag_i_dn*c_j_dn_cdag_j_up)+(c_i_up_cdag_j_up*cdag_i_dn_c_j_dn)

                    sx3=(c_i_dn_cdag_i_up*c_j_up_cdag_j_dn)+(c_i_dn_cdag_j_up*cdag_i_up_c_j_up)

                    sx4=(c_i_up_cdag_i_dn*c_j_up_cdag_j_dn)+(c_i_up_cdag_j_dn*cdag_i_dn_c_j_up)
                   
           
                    sx_i_sx_j=sx_i_sx_j+(sx1+sx2+sx3+sx4)*(1/4.0d0)*((-1)**(abs(i-j)))
                    sy_i_sy_j=sy_i_sy_j+(-sx1+sx2+sx3-sx4)*(1/4.0d0)*((-1)**(abs(i-j)))

                endif
            enddo
        endif
    enddo

    !print*,'Sx_Sx_a',sx_i_sx_j !/(nos/3)
    !print*,'Sy_Sy_a',sy_i_sy_j !/(nos/3)
endsubroutine s_xy_corr_a

subroutine sz_corr_a(temp2)
    use array
    use global
    implicit none
    integer::j1,l
    doubleprecision::c_i_up_cdag_i_up,c_j_up_cdag_j_up,c_i_dn_cdag_j_up,c_j_dn_cdag_j_dn,cdag_i_dn_c_j_up 
    double precision::cdag_i_up_c_j_up,c_i_dn_cdag_j_dn,cdag_i_up_c_j_dn,c_i_up_cdag_j_dn,c_i_dn_cdag_i_dn
    doubleprecision::c_i_up_cdag_j_up,fermi_fn,sz1,sz2,sz3,sz4,cdag_i_dn_c_j_dn,temp2

    sz_i_sz_j=0.0d0
            
    do i=1,nos
        if(tag(i).eq.1)then
            do j=1,nos
                if((tag(j).eq.1).and.(i.eq.j))then
    
                    c_i_up_cdag_i_up=0.0d0
                    c_j_up_cdag_j_up=0.0d0
                    c_i_up_cdag_j_up=0.0d0
                    c_i_dn_cdag_i_dn=0.0d0
                    c_j_dn_cdag_j_dn=0.0d0
                    c_i_dn_cdag_j_up=0.0d0
                    c_i_up_cdag_j_dn=0.0d0
                    c_i_dn_cdag_j_dn=0.0d0
                    cdag_i_up_c_j_up=0.0d0
                    cdag_i_dn_c_j_up=0.0d0
                    cdag_i_up_c_j_dn=0.0d0
                    cdag_i_dn_c_j_dn=0.0d0   

                    do l=1,2*nos

                        fermi_fn=(1/(1.0d0+(exp((evl_s(l)-mu)/temp2))))

                        c_i_up_cdag_i_up=c_i_up_cdag_i_up+((h(2*i-1,l))*conjg(h(2*i-1,l)))*(1-fermi_fn)   !1

                        c_j_up_cdag_j_up=c_j_up_cdag_j_up+((h(2*j-1,l))*conjg(h(2*j-1,l)))*(1-fermi_fn)    !2

                        c_i_up_cdag_j_up=c_i_up_cdag_j_up+((h(2*i-1,l))*conjg(h(2*j-1,l)))*(1-fermi_fn)    !3

                        c_i_dn_cdag_i_dn=c_i_dn_cdag_i_dn+((h(2*i,l))*conjg(h(2*i,l)))*(1-fermi_fn)    !4

                        c_j_dn_cdag_j_dn=c_j_dn_cdag_j_dn+((h(2*j-1,l))*conjg(h(2*j-1,l)))*(1-fermi_fn)   !5

                        c_i_dn_cdag_j_up=c_i_dn_cdag_j_up+((h(2*i,l))*conjg(h(2*j-1,l)))*(1-fermi_fn)     !6

                        c_i_up_cdag_j_dn=c_i_up_cdag_j_dn+((h(2*i-1,l))*conjg(h(2*j,l)))*(1-fermi_fn)     !7

                        c_i_dn_cdag_j_dn=c_i_dn_cdag_j_dn+((h(2*i,l))*conjg(h(2*j,l)))*(1-fermi_fn)      !8

                        cdag_i_up_c_j_up=cdag_i_up_c_j_up+((h(2*j-1,l))*conjg(h(2*i-1,l)))*fermi_fn       !9

                        cdag_i_dn_c_j_up=cdag_i_dn_c_j_up+((h(2*j-1,l))*conjg(h(2*i,l)))*fermi_fn   !10

                        cdag_i_up_c_j_dn=cdag_i_up_c_j_dn+((h(2*j,l))*conjg(h(2*i-1,l)))*fermi_fn     !11

                        cdag_i_dn_c_j_dn=cdag_i_dn_c_j_dn+((h(2*j-1,l))*conjg(h(2*i,l)))*fermi_fn   !12
                    enddo
        
                    sz1=(c_i_up_cdag_i_up*c_j_up_cdag_j_up)+(c_i_up_cdag_j_up*cdag_i_up_c_j_up)
                    
                    sz2=(c_i_dn_cdag_i_dn*c_j_up_cdag_j_up)+(c_i_dn_cdag_j_up*cdag_i_dn_c_j_up)
                    
                    sz3=(c_i_up_cdag_i_up*c_j_dn_cdag_j_dn)+(c_i_up_cdag_j_dn*cdag_i_up_c_j_dn)
                    
                    sz4=(c_i_dn_cdag_i_dn*c_j_dn_cdag_j_dn)+(c_i_dn_cdag_j_dn*cdag_i_dn_c_j_dn)

                    
           
                    sz_i_sz_j=sz_i_sz_j+(sz1-sz2-sz3+sz4)*((-1)**(abs(i-j)))*(1/4.0d0)
        
                endif
            enddo
        endif
    enddo

endsubroutine sz_corr_a



subroutine bloch_states(eignv_mat)
    use global
    use array
    implicit none 
    integer::k,level,comp,p,m
    double precision::k_r,exp_k_r,norm,ky,kx,pi
    double precision::bottom_band((d/2)**2),flat_band((d/2)**2),top_band((d/2)**2),dummy_value,eps_
    complex*16::eignv_mat(2*nos,2*nos)

    allocate(mat_1(((d/2)**2),2*nos),mat_2(((d/2)**2),2*nos),mat_3(((d/2)**2),2*nos),vec(2*nos))
    allocate(mat_1_up(((d/2)**2),2*nos),mat_2_up(((d/2)**2),2*nos),mat_3_up(((d/2)**2),2*nos))
    allocate(mat_1_down(((d/2)**2),2*nos),mat_2_down(((d/2)**2),2*nos),mat_3_down(((d/2)**2),2*nos))
    
    eps_=10!1e-6
    pi=4.0*atan(1.0)
    dummy_value=-123456789.0
    bottom_band=dummy_value
    flat_band=dummy_value
    top_band=dummy_value  

    mat_1=complex(0.0d0,0.0d0)
    mat_2=complex(0.0d0,0.0d0)
    mat_3=complex(0.0d0,0.0d0) 
    mat_1_up=complex(0.0d0,0.0d0) 
    mat_2_up=complex(0.0d0,0.0d0) 
    mat_3_up=complex(0.0d0,0.0d0) 
    mat_1_down=complex(0.0d0,0.0d0) 
    mat_2_down=complex(0.0d0,0.0d0) 
    mat_3_down=complex(0.0d0,0.0d0) 
    
    do m=1,2*nos
        e=evl_s(m)
        level=m
        i=0
        do k_y=1,(d/2)  
            ky =(2*pi/(d/2))*k_y
            do k_x=1,(d/2)
                kx=(2*pi/(d/2))*k_x
                k=k_x+(k_y-1)*(d/2)
                i=i+1    
                do ly=0,(d/2)-1  
                    do lx=0,(d/2)-1 
                        do a1=0,2
                            do sigma=0,1           
                                comp=(lx)*6+(a1)*2+sigma+(ly)*6*(d/2)              
                                vec(comp+1)=complex(0.0d0,0.0d0)    
                                do jy=0,(d/2)-1  
                                    do jx=0,(d/2)-1 
                                        nx=(mod((lx+jx),(d/2)))
                                        ny=(mod((ly+jy),(d/2)))                
                                        i_3=(nx)*6+(a1)*2+sigma+(ny)*6*(d/2)          
                                        k_r=(kx*(jx)+ky*(jy))
                                        exp_k_r=complex(cos(k_r),-sin(k_r)) 
                                        vec(comp+1)=vec(comp+1)+exp_k_r*(eignv_mat( (i_3+1), m)) 
                                    enddo          
                                enddo
                            enddo      
                        enddo        
                    enddo          
                enddo    
      
                norm=0.0d0 
                do p=1,2*nos
                    norm=norm+vec(p)*conjg(vec(p)) 
                enddo  ! p loop 
               
             
                if(norm.gt.eps_)then            
                    if((bottom_band(k) .eq.dummy_value))then 
                        bottom_band(k)=e
                        do p=1,2*nos
                            mat_1(k,p)=vec(p)  
                            if(mod(p,2).ne.0) mat_1_up(k,p)=vec(p)
                            if(mod(p,2).eq.0) mat_1_down(k,p)=vec(p)                
                        enddo
                    elseif((flat_band(k) .eq.dummy_value).and. (abs(e-bottom_band(k)).gt.eps_))then
                        flat_band(k)=e
                        do p=1,2*nos               
                            mat_2(k,p)=vec(p)                
                            if(mod(p,2).ne.0) mat_2_up(k,p)=vec(p)
                            if(mod(p,2).eq.0) mat_2_down(k,p)=vec(p) 
                        enddo              
                    ! elseif((top_band(k) .eq.dummy_value).and.(abs(e-flat_band(k)).gt.eps_).and.(abs(e-bottom_band(k)).gt.eps_))then
                    !     top_band(k)=e
                    !     do p=1,2*nos
                    !         mat_3(k,p)=vec(p)  
                    !         if(mod(p,2).ne.0) mat_3_up(k,p)=vec(p)
                    !         if(mod(p,2).eq.0) mat_3_down(k,p)=vec(p)
                    !     enddo             
                    endif        
                endif              
            enddo           
        enddo         
    enddo  
      
    do k_x=1,d/2
        kx =(2*pi/(d/2))*k_x
        do k_y=1,d/2  
            ky=(2*pi/(d/2))*k_y
            k=k_x+(k_y-1)*(d/2)
            write(101,*)kx,ky,k,bottom_band(k),flat_band(k),top_band(k)     
        enddo
    enddo
end subroutine bloch_states

subroutine bloch_states6(eignv_mat)
    use global
    use array
    implicit none 
    integer::k,level,comp,p,m
    double precision::k_r,exp_k_r,norm,ky,kx,pi
    double precision::bottom_band1((d/2)**2),flat_band1((d/2)**2),top_band1((d/2)**2),dummy_value,eps_
    double precision::bottom_band2((d/2)**2),flat_band2((d/2)**2),top_band2((d/2)**2)
    complex*16::eignv_mat(2*nos,2*nos)
    allocate(mat_1(((d/2)**2),2*nos),mat_2(((d/2)**2),2*nos),mat_3(((d/2)**2),2*nos),vec(2*nos))
    allocate(mat_4(((d/2)**2),2*nos),mat_5(((d/2)**2),2*nos),mat_6(((d/2)**2),2*nos))
  
    allocate(mat_1_up(((d/2)**2),2*nos),mat_2_up(((d/2)**2),2*nos),mat_3_up(((d/2)**2),2*nos))
    allocate(mat_4_up(((d/2)**2),2*nos),mat_5_up(((d/2)**2),2*nos),mat_6_up(((d/2)**2),2*nos))
  
    allocate(mat_1_down(((d/2)**2),2*nos),mat_2_down(((d/2)**2),2*nos),mat_3_down(((d/2)**2),2*nos))
    allocate(mat_4_down(((d/2)**2),2*nos),mat_5_down(((d/2)**2),2*nos),mat_6_down(((d/2)**2),2*nos))
  
    !vec=complex(0.0d0,0.0d0)
    pi=4.0*atan(1.0)
    eps_=1e-6
    dummy_value=-123456789.0
    bottom_band1=dummy_value
    bottom_band2=dummy_value
    flat_band1=dummy_value
    flat_band2=dummy_value
    top_band1=dummy_value  
    top_band2=dummy_value  
  
    ! print*,'here 1' 
    mat_1=complex(0.0d0,0.0d0)
    mat_2=complex(0.0d0,0.0d0)
    mat_3=complex(0.0d0,0.0d0) 
    mat_4=complex(0.0d0,0.0d0)
    mat_5=complex(0.0d0,0.0d0)
    mat_6=complex(0.0d0,0.0d0) 
    mat_1_up=complex(0.0d0,0.0d0) 
    mat_2_up=complex(0.0d0,0.0d0) 
    mat_3_up=complex(0.0d0,0.0d0) 
    mat_4_up=complex(0.0d0,0.0d0) 
    mat_5_up=complex(0.0d0,0.0d0) 
    mat_6_up=complex(0.0d0,0.0d0) 
    mat_1_down=complex(0.0d0,0.0d0) 
    mat_2_down=complex(0.0d0,0.0d0) 
    mat_3_down=complex(0.0d0,0.0d0) 
    mat_4_down=complex(0.0d0,0.0d0) 
    mat_5_down=complex(0.0d0,0.0d0) 
    mat_6_down=complex(0.0d0,0.0d0) 
  
    do m=1,2*nos
      e=evl_s(m)
      level=m
      i=0
      do k_y=1,(d/2)  !d2
       !  print*,'here 2'  
       ky =(2*pi/(d/2))*k_y
       do k_x=1,(d/2)
         kx=(2*pi/(d/2))*k_x
         ! print*,'for',p1,k1 
         k=k_x+(k_y-1)*(d/2)	 
         i=i+1   
  
         do ly=0,(d/2)-1  
            do lx=0,(d/2)-1 
                do a1=0,2
                    do sigma=0,1           
                        comp=(lx)*6+(a1)*2+sigma+(ly)*6*(d/2)              
                        vec(comp+1)=complex(0.0d0,0.0d0)    
                        do jy=0,(d/2)-1  
                            do jx=0,(d/2)-1 
                                nx=(mod((lx+jx),(d/2)))
                                ny=(mod((ly+jy),(d/2)))                
                                i_3=(nx)*6+(a1)*2+sigma+(ny)*6*(d/2)          
                                k_r=(kx*(jx)+ky*(jy))
                                exp_k_r=complex(cos(k_r),-sin(k_r)) 
                                vec(comp+1)=vec(comp+1)+exp_k_r*(eignv_mat( (i_3+1), m)) 
                            enddo          
                        enddo
                    enddo      
                enddo        
            enddo          
        enddo  
        
  
          
  
    
         norm=0.0d0 
         do p=1,2*nos
            norm=norm+vec(p)*conjg(vec(p)) 
         enddo  ! p loop 
               
          !print*,m,e,k_x,k_y,norm             
          if(norm.gt.eps_)then
            
            if((bottom_band1(k) .eq.dummy_value))then 
              bottom_band1(k)=e
             
              do p=1,2*nos
                mat_1(k,p)=vec(p)
               
                
                if(mod(p,2).ne.0) mat_1_up(k,p)=vec(p)
                if(mod(p,2).eq.0) mat_1_down(k,p)=vec(p)
  
                
              enddo
              
  
            elseif((bottom_band2(k) .eq.dummy_value).and. (abs(e-bottom_band1(k)).gt.eps_))then
              bottom_band2(k)=e
              !print*,k,e,m,'flat_band'
              do p=1,2*nos
                
                mat_2(k,p)=vec(p)
            
                if(mod(p,2).ne.0) mat_2_up(k,p)=vec(p)
                if(mod(p,2).eq.0) mat_2_down(k,p)=vec(p)
  
              enddo
  
             ! print*,k,e,kx,ky
            elseif((flat_band1(k) .eq.dummy_value).and. (abs(e-bottom_band1(k)).gt.eps_)&
              .and. (abs(e-bottom_band2(k)).gt.eps_))then
              flat_band1(k)=e
              !print*,k,e,m,'flat_band'
              do p=1,2*nos
                
                mat_3(k,p)=vec(p)
                
                if(mod(p,2).ne.0) mat_3_up(k,p)=vec(p)
                if(mod(p,2).eq.0) mat_3_down(k,p)=vec(p)
  
              enddo
              
            elseif((flat_band2(k) .eq.dummy_value).and. (abs(e-bottom_band1(k)).gt.eps_)&
              .and. (abs(e-bottom_band2(k)).gt.eps_).and. (abs(e-flat_band1(k)).gt.eps_))then
              flat_band2(k)=e
              !print*,k,e,m,'flat_band'
              do p=1,2*nos
                
                mat_4(k,p)=vec(p)
                
                if(mod(p,2).ne.0) mat_4_up(k,p)=vec(p)
                if(mod(p,2).eq.0) mat_4_down(k,p)=vec(p)
  
              enddo
              
            elseif((top_band1(k) .eq.dummy_value).and. (abs(e-bottom_band1(k)).gt.eps_)&
              .and. (abs(e-bottom_band2(k)).gt.eps_).and. (abs(e-flat_band1(k)).gt.eps_)&
              .and. (abs(e-flat_band2(k)).gt.eps_))then
              top_band1(k)=e
              !print*,k,e,m,'flat_band'
              do p=1,2*nos
                
                mat_5(k,p)=vec(p)
                write(80,*)k,p,mat_5(k,p)
                if(mod(p,2).ne.0) mat_5_up(k,p)=vec(p)
                if(mod(p,2).eq.0) mat_5_down(k,p)=vec(p)
  
              enddo
  
            elseif((top_band2(k) .eq.dummy_value).and. (abs(e-bottom_band1(k)).gt.eps_)&
              .and. (abs(e-bottom_band2(k)).gt.eps_).and. (abs(e-flat_band1(k)).gt.eps_)&
              .and. (abs(e-flat_band2(k)).gt.eps_).and. (abs(e-top_band1(k)).gt.eps_))then
              top_band2(k)=e
              !print*,k,e,m,'flat_band'
              do p=1,2*nos
                
                mat_6(k,p)=vec(p)
                write(81,*)k,p,mat_6(k,p)
                if(mod(p,2).ne.0) mat_6_up(k,p)=vec(p)
                if(mod(p,2).eq.0) mat_6_down(k,p)=vec(p)
  
              enddo
      
            endif        
         endif              
        enddo           
      enddo         
    enddo  
  
    
    do k_x=1,d/2
        kx =(2*pi/(d/2))*k_x
        do k_y=1,d/2  
            ky=(2*pi/(d/2))*k_y
            k=k_x+(k_y-1)*(d/2)
            write(101,*)kx,ky,k,bottom_band1(k),bottom_band2(k),flat_band1(k),flat_band2(k),top_band1(k),top_band2(k)      
            flush(101)
        enddo
    enddo
    
   
end subroutine bloch_states6

subroutine bottindex(i1,i2,trace)
    use global 
    use array
    use s_vd
    implicit none
    integer::l,k,n,i1,nfill,i2
    doubleprecision::pi,l1,k1,x2,y2,z2,phi1,theta1,det,det1,trace,mu1,mu2
    complex*16::proj,proj_op(2*nos,2*nos),check(2*nos,2*nos),e_theta(2*nos,2*nos),e_phi(2*nos,2*nos),u1(2*nos,2*nos),w1(2*nos,2*nos)
    complex*16::u2(2*nos,2*nos),w2(2*nos,2*nos),prod(2*nos,2*nos),u3(2*nos,2*nos),w3(2*nos,2*nos)
    complex*16::u4(2,2),FindDet,proj_op1(2*nos,2*nos),det2,prod1(2*nos,2*nos),prod2(2*nos,2*nos)

    pi=4.0*atan(1.0)
    proj_op=complex(0.0d0,0.0d0)
    mu1=0.410d0
    mu2=2.40d0
    nfill=2*nos/6
    ! nfill=0
    ! do i=1,2*nos
    !    if (evl_s(i).le.mu) nfill= nfill+1
    ! enddo
    ! print*,nfill

    do i=1,2*nos    
        do j=1,2*nos 
            proj=cmplx(0.0d0,0.0d0)          
            do l=(i1-1)*nfill+1,(i2-1)*nfill
              
                proj=proj+h(i,l)*conjg(h(j,l))
                
            enddo
            proj_op(i,j)=proj
        enddo
    enddo
     

    do i=1,nos           
        if(tag(i).eq.1)then
            phi(i)=2*pi*x(i)/d
            theta(i)=2*pi*y(i)/d
        endif

        if(tag(i).eq.2)then 
            phi(i)=2*pi*x(i)/d
            theta(i)=2*pi*y(i)/d
        endif      

        if(tag(i).eq.3)then
            phi(i)=2*pi*x(i)/d
            theta(i)=2*pi*y(i)/d 
        endif     
    enddo


    do i=1,nos
        x2=(2+cos(theta(i)))*cos(phi(i))
        y2=(2+cos(theta(i)))*sin(phi(i))
        z2=sin(theta(i))
        write(57,*) phi(i),theta(i),x2,y2,z2
    enddo


    do i=1,nos
        e_theta(2*i-1,2*i-1)=exp(cmplx(0.0d0,theta(i)))
        e_theta(2*i,2*i)=exp(cmplx(0.0d0,theta(i)))
        e_phi(2*i-1,2*i-1)=exp(cmplx(0.0d0,phi(i)))
        e_phi(2*i,2*i)=exp(cmplx(0.0d0,phi(i)))
    enddo
 
    ! do i=1,2*nos
    !     write(30,*) ((e_theta(i,j)),j=1,2*nos)
    !     write(31,*) ((e_phi(i,j)),j=1,2*nos)
    ! enddo

    u3=matmul(proj_op,e_theta)
    u1=matmul(u3,proj_op)

    w3=matmul(proj_op,e_phi)
    w1=matmul(w3,proj_op)
 
    ! do i=1,2*nos
    !     write(34,*) ((u1(i,j)),j=1,2*nos)
    !     write(35,*) ((w1(i,j)),j=1,2*nos)
    ! enddo

    do i=1,2*nos
        do j=1,2*nos
            u2(i,j)=conjg(u1(j,i)) 
            w2(i,j)=conjg(w1(j,i)) 
        enddo
    enddo
 

    ! u2=transpose(conjg(u1))
    ! w2=transpose(conjg(w1))

    prod1=matmul(u1,w1)
    prod2=matmul(prod1,u2)
    prod=matmul(prod2,w2)

    
    ! call s_v_d(u1,2*nos,2*nos)

    ! do i=1,2*nos
    !     do j=1,2*nos
    !         write(3000,*) u_1(i,j),vt(i,j)  
    !     enddo
    ! enddo    

    call diagonalization1(prod)



    trace=0.0d0

    do i=1,2*nos
        if(abs(real(alpha(i))).gt.1e-6)then
            trace=trace+aimag(log(alpha(i)))/(2*pi)
        endif
       !print*,alpha(i),log(alpha(i))
    enddo

    
    ! u4(1,1)=cmplx(2,3)
    ! u4(1,2)=cmplx(0,5)
    ! u4(2,1)=cmplx(0,1)
    ! u4(2,2)=3

    ! u4=transpose(conjg(u4))
    !print*,u4(1,2),u4(2,1)
    ! det1=finddet(prod,2*nos)
    ! print*,det1

    
endsubroutine bottindex



subroutine struc_factor
    use array
    use global
    implicit none
    double precision::pi,m_x,m_y,m_z,mimj,qx,qy,rx,ry,sf_a,sf_b_c
    pi=4.0*atan(1.0)
    qx=0.0d0
    qy=0.0d0
    !s_f=0.0d0
    sf_A=0.0d0
    
    do i=1,nos ! nos is no of sites  
        if(tag(i).eq.1)then
            do j=1,nos
                if(tag(j).eq.1)then
                    
                    m_x=(mi(i)*sin(th(i))*cos(ph(i)))*(mi(j)*sin(th(j))*cos(ph(j)))
                    m_y=(mi(i)*sin(th(i))*sin(ph(i)))*(mi(j)*sin(th(j))*sin(ph(j)))
                    m_z=(mi(i)*cos(ph(i)))*(mi(j)*cos(ph(j)))
                    mimj=m_x+m_y+m_z
                    rx=qx*(x(i)-x(j))
                    ry=qy*(y(i)-y(j))
       
                    sf_A=sf_A+exp(cmplx(0.0d0,rx+ry))*(mimj)        

                endif
            enddo
        endif
    enddo
    sf_B_C=0.0d0
    do i=1,nos 
        if((tag(i).eq.2).or.(tag(i).eq.3))then
            do j=1,nos
                if((tag(j).eq.2).or.(tag(j).eq.3))then                  
                    m_x=(mi(i)*sin(th(i))*cos(ph(i)))*(mi(j)*sin(th(j))*cos(ph(j)))
                    m_y=(mi(i)*sin(th(i))*sin(ph(i)))*(mi(j)*sin(th(j))*sin(ph(j)))
                    m_z=(mi(i)*cos(ph(i)))*(mi(j)*cos(ph(j)))
                    mimj=m_x+m_y+m_z
                    rx=qx*(x(i)-x(j))
                    ry=qy*(y(i)-y(j))
                    sf_B_C=sf_B_C+exp(cmplx(0.0d0,rx+ry))*(mimj)        

                endif
            enddo
        endif
    enddo
endsubroutine struc_factor


subroutine proj_dos(temp2)
    use global
    use array
    implicit none
    integer::nwp
    double precision::eta_,wmin,wmax,dw,w,p_dos_a,fermi_fn,p_dos_b,p_dos_c,mat_a,mat_b,mat_c,pi,temp2
    
    pi=4.0*atan(1.0) 
    wmin=-20
    wmax=20
    nwp=2000
    dw=abs(wmax-wmin)/nwp
    w=wmin
    eta_=0.050d0

    do ie=1,nwp
        
        w=w+dw
        p_dos_a=0.0d0
        p_dos_b=0.0d0
        p_dos_c=0.0d0
        
        
        
        do i=1,nos
            if(tag(i).eq.1)then
                do j=1,2*nos
                    fermi_fn=(1/(1.0d0+(exp((evl_s(j)-mu)/temp2))))
                    mat_a=((h(2*i,j)*conjg(h(2*i,j)))+(h(2*i-1,j)*conjg(h(2*i-1,j))))
                    p_dos_a=p_dos_a+((eta_/pi)/((w-(evl_s(j)-mu))**2+((eta_)**2)))*mat_a
                enddo
            endif
        enddo
        


        do i=1,nos       
            if (tag(i).eq.2)then
                do j=1,2*nos
                    fermi_fn=(1/(1.0d0+(exp((evl_s(j)-mu)/temp2))))
                    mat_b=((h(2*i,j)*conjg(h(2*i,j)))+(h(2*i-1,j)*conjg(h(2*i-1,j))))
                    p_dos_b=p_dos_b+((eta_/pi)/((w-(evl_s(j)-mu))**2+((eta_)**2)))*mat_b
                enddo
            endif
        enddo



        do i=1,nos
            if (tag(i).eq.3)then
                do j=1,2*nos
                    fermi_fn=(1/(1.0d0+(exp((evl_s(j)-mu)/temp2))))
                    mat_c=((h(2*i,j)*conjg(h(2*i,j)))+(h(2*i-1,j)*conjg(h(2*i-1,j))))
                    p_dos_c=p_dos_c+((eta_/pi)/((w-(evl_s(j)-mu))**2+((eta_)**2)))*mat_c
                enddo
            endif
        enddo
        write(300,*)w,p_dos_a/(2*nos),p_dos_b/(2*nos),p_dos_c/(2*nos),(p_dos_a+p_dos_b+p_dos_c)/(2*nos)       
                                    
    enddo 
endsubroutine proj_dos


subroutine spec_func
    use array
    use global
    implicit none
    double precision::mat_e,c1,pi,eta,spec_f,kx,ky
    complex*16::c2
    complex*16::mat_sum(2*nos,2*nos)
    integer::l
    pi=4.0*atan(1.0)
    eta=0.025000d0   
    e=2.0d0  
    mat_sum=0.0d0
    do l=1,nos     
        do j=1,nos
            do b1=1,2*nos
                mat_e=((h(2*j,b1)*conjg(h(2*l,b1)))+(h(2*j-1,b1)*conjg(h(2*l-1,b1))))
                mat_sum(l,j)=mat_sum(l,j)+mat_e*((eta/(2*pi))/((e-evl_s(b1))**2+((eta/2)**2)))                     
            enddo
        enddo
    enddo
          
    do k_y=1,(d/2)   
        ky =(2*pi/(d/2))*k_y
        do k_x=1,(d/2)
            kx=(2*pi/(d/2))*k_x
            spec_f=0.0d0 
            do l=1,nos     
                do j=1,nos                                
                    c1=(kx*(x(j)-x(l))+ky*(y(j)-y(l)))
                    c2=complex(cos(c1),sin(c1))
                    spec_f=spec_f+c2*mat_sum(l,j)                      
                enddo   
            enddo    
            write(120,*) kx,ky,spec_f/nos
            flush(120)
        enddo       
    enddo
end subroutine spec_func


subroutine occupation(temp2)
    use array
    use global
    implicit none
    integer::l,k,a,b
    doubleprecision::fermi_fn,temp2,n_a_up,n_b_up,n_c_up,n_a_dn,n_b_dn,n_c_dn

    n_a_up=0.0d0
    n_b_up=0.0d0
    n_c_up=0.0d0
    n_a_dn=0.0d0
    n_b_dn=0.0d0
    n_c_dn=0.0d0

    do l=1,nos
        if(tag(l).eq.1)then
            a=2*l-1
            do k=1,2*nos
                fermi_fn=(1/(1.0d0+(exp((evl_s(k)-mu)/temp2))))         
                b=k
                n_a_up=n_a_up+(h(a,b)*conjg(h(a,b)))*fermi_fn
                n_a_dn=n_a_dn+(h(a+1,b)*conjg(h(a+1,b)))*fermi_fn
            enddo
        endif
          
        if((tag(l).eq.2))then
            a=2*l-1
            do k=1,2*nos
                fermi_fn=(1/(1.0d0+(exp((evl_s(k)-mu)/temp2))))
                b=k
                n_b_up=n_b_up+(h(a,b)*conjg(h(a,b)))*fermi_fn
                n_b_dn=n_b_dn+(h(a+1,b)*conjg(h(a+1,b)))*fermi_fn
            enddo
        endif

        if(tag(l).eq.3)then
            a=2*l-1
            do k=1,2*nos
                fermi_fn=(1/(1.0d0+(exp((evl_s(k)-mu)/temp2))))
                b=k
                n_c_up=n_c_up+(h(a,b)*conjg(h(a,b)))*fermi_fn
                n_c_dn=n_c_dn+(h(a+1,b)*conjg(h(a+1,b)))*fermi_fn
            enddo
        endif
    enddo
    write(800,*)n_a_up/unit_cells,n_a_dn/unit_cells,n_b_up/unit_cells,n_b_dn/unit_cells,n_c_up/unit_cells,n_c_dn/unit_cells
endsubroutine occupation


subroutine chern_number
    use global 
    use array
    implicit none
    integer::p,count_k,kp1,p_,kyp1,kxp1,kp2,kp1p2
    complex*16::f12,u1,u2,u3,u4,u1k,u2k,u1k2,u2k1
    complex*16::nk_up((2*nos)),nkp1_up((2*nos)),nkp2_up((2*nos)),nkp1p2_up((2*nos))
    double precision::k_r,kp1_r,kp2_r,kp1p2_r,norm_k,norm_kp1,norm_kp2,norm_kp1p2
    !complex*16::nkp1_down((2*nos)),nkp2_down((2*nos)),nkp1p2_down((2*nos)),nk_down((2*nos))
    double precision::pi,ky,kx

    f12=complex(0.0d0,0.0d0)
     
    count_k=1
    pi=4.0*atan(1.0) 
    do k_y=1,d/2
        ky =(2*pi/(d/2))*k_y
        do k_x=1,d/2
            kx=(2*pi/(d/2))*k_x
       
            kxp1=(mod(k_x,d/2))+1
            kyp1=(mod(k_y,d/2))+1
  
            kp1=kxp1+(k_y-1)*(d/2)    ! along kx direction
            kp2=k_x+(kyp1-1)*(d/2)    ! along ky direction
            kp1p2=kxp1+(kyp1-1)*(d/2) 
                 
            norm_k=0.0d0
            norm_kp1=0.0d0
            norm_kp2=0.0d0
            norm_kp1p2=0.0d0       
  
            do p=1,2*nos
  
                norm_k=norm_k+(mat_temp(count_k,p)*conjg(mat_temp(count_k,p)))
                norm_kp1=norm_kp1+(mat_temp(kp1,p)*conjg(mat_temp(kp1,p)))
                norm_kp2=norm_kp2+(mat_temp(kp2,p)*conjg(mat_temp(kp2,p)))
                norm_kp1p2=norm_kp1p2+(mat_temp(kp1p2,p)*conjg(mat_temp(kp1p2,p)))
  
            enddo
  
            do ly=1,d/2  
                do lx=1,d/2 
  
                    k_r=kx*lx+ky*ly 
                    kp1_r=((2*pi/(d/2))*kxp1)*lx + ((2*pi/(d/2))*k_y)*ly
                    kp2_r=((2*pi/(d/2))*k_x)*lx + ((2*pi/(d/2))*kyp1)*ly
                    kp1p2_r=((2*pi/(d/2))*kxp1)*lx + ((2*pi/(d/2))*kyp1)*ly
  
  
                    do a1=1,3
                        do sigma=1,2
                            ! if(a1.le.2)  p_=(lx-1)*4+(a1-1)*2+sigma+(ly-1)*6*(d/2) 
                            ! if(a1.eq.3)  p_=4*(d/2)+sigma+2*(lx-1)+(ly-1)*6*(d/2)   
                            p_=(lx-1)*6+(a1-1)*2+sigma+(ly-1)*6*(d/2)

                            nk_up(p_)=complex(cos(k_r),sin(k_r))*(mat_temp(count_k,p_)/sqrt(norm_k))
                            nkp1_up(p_)=complex(cos(kp1_r),sin(kp1_r))*(mat_temp(kp1,p_)/sqrt(norm_kp1))
                            nkp2_up(p_)=complex(cos(kp2_r),sin(kp2_r))*(mat_temp(kp2,p_)/sqrt(norm_kp2))
                            nkp1p2_up(p_)=complex(cos(kp1p2_r),sin(kp1p2_r))*(mat_temp(kp1p2,p_)/sqrt(norm_kp1p2)) 
                            ! print*,p_,nkp1_up(p_)
                        enddo
                    enddo
                enddo
            enddo 
              
            u1=complex(0.0d0,0.0d0) 
            u2=complex(0.0d0,0.0d0) 
            u3=complex(0.0d0,0.0d0) 
            u4=complex(0.0d0,0.0d0) 
  
            do p=1,(2*nos) 
                u1=u1+((conjg(nk_up(p))*nkp1_up(p)))     
                u2=u2+((conjg(nkp1_up(p))*nkp1p2_up(p)))       
                u3=u3+((conjg(nkp2_up(p))*nkp1p2_up(p)))      
                u4=u4+((conjg(nk_up(p))*nkp2_up(p)))
            enddo
  
            u1k=u1/sqrt(u1*conjg(u1))
            u2k1=u2/sqrt(u2*conjg(u2))
            u1k2=u3/sqrt(u3*conjg(u3))
            u2k=u4/sqrt(u4*conjg(u4))
  
  
            f12=f12+log(u1k*u2k1*(1/u1k2)*(1/u2k))
       
            count_k=count_k+1
            !write(110,*) kx,ky,real(f12/(2*pi*complex(0.0d0,1.0d0)))
  
         enddo
    enddo
    
    ch_no=aimag(f12)/(2*pi)
    
    print*,ch_no
    
endsubroutine chern_number

complex*16 function FindDet(matrix, n)
    implicit none
    complex*16,dimension(n,n):: matrix
    integer,intent(in) :: n
    complex*16::m,temp
    integer::i,j,k,l
    logical:: DetExists = .true.
    l = 1
    !Convert to upper triangular form
    do k = 1, n-1
        if (matrix(k,k) == 0)then          
            DetExists = .false.
            do i = k+1, n
                if (matrix(i,k) /= 0)then
                    do j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    enddo
                    DetExists = .true.
                    l=-l
                    EXIT
                endif
            enddo
            if (DetExists .eqv. .false.)then
               
                FindDet = 0
                return
            endif
        endif
        do j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            do i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            enddo
        enddo
    enddo
   
    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    do i = 1, n
        FindDet = FindDet * matrix(i,i)
    enddo
   
endfunction FindDet

function zdet(A) result(x)
    ! compute the determinant of a real matrix using an LU factorization
    complex*16, intent(in) :: A(:, :)
    complex(16) :: x
    integer :: i
    ! LAPACK variables:
    integer :: info, n
    integer, allocatable :: ipiv(:)
    complex*16, allocatable :: At(:,:)

    n = size(A(1,:))
    !call assert_shape(A, [n, n], "det", "A")
    allocate(At(n,n), ipiv(n))
    At = A
    call zgetrf(n, n, At, n, ipiv, info)
    if(info /= 0) then
       print *, "zgetrf returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, and division by zero will occur if it is used"
          print *, "to solve a system of equations."
       end if
       !call stop_error('inv: zgetrf error')
    end if

    ! for details on the computation, compare the comment in ddet().
    x = 1.0d0 + 0*i_
    do i = 1,n
       if(ipiv(i) /= i) then  ! additional sign change
          x = -x*At(i,i)
       else
          x = x*At(i,i)
       endif
    end do
  end function zdet


  complex*16 function DET2(aa,n2)
  use array
  use global
  implicit none
  complex*16::aa(n2,n2)
  complex*16::tmp,c(size(aa,dim=1),size(aa,dim=2))
  complex*16::max
  integer::n2,n,k,l,m,num(size(aa,dim=1))
  
      n=size(aa,dim=1)
      det2=1
      do k=1,n
          max=aa(k,k);num(k)=k;
          do i=k+1,n 
              if(abs(max)<abs(aa(i,k))) then
                  max=aa(i,k)
                  num(k)=i
              endif
          enddo
          if (num(k)/=k) then
              do l=k,n 
                  tmp=aa(k,l)
                  aa(k,l)=aa(num(k),l)
                  aa(num(k),l)=tmp
              enddo
              det2=-1*det2
          endif
          do m=k+1,n
              c(m,k)=aa(m,k)/aa(k,k)
              do l=k,n 
                  aa(m,l)=aa(m,l)-c(m,k)*aa(k,l)
              enddo
          enddo !There we made matrix triangular!	
      enddo
  
      do i=1,n
      det2=det2*aa(i,i)
      enddo
      return
end function

!--------------------------------------------
subroutine s_v_d(d_eff,m1_,n1)
    use s_vd
    use global
    use array
    implicit none
    INTEGER::m1_,n1
    iNTEGER::LDA, LDU, LDVT
    INTEGER::LWMAX,lwork
    INTEGER::INFO
    complex*16::d_eff(2*nos,2*nos)
    
    LDA = M1_
    LDU = M1_
    LDVT = N1 
    LWMAX = 3*m1_

    allocate( U_1( LDU, M1_ ),VT( LDVT, N1 ),WORK( LWMAX ) ,sigma1( M1_ ),RWORK( 5*M1_ ))

    LWORK = -1
    CALL ZGESVD( 'All', 'All', m1_, N1,d_eff, LDA, sigma1, U_1, LDU, VT, LDVT,WORK, LWORK, RWORK, INFO )
    
    LWORK = MIN( LWMAX,INT(WORK(1)))
    !*     Compute SVD.

    CALL ZGESVD( 'All', 'All', M1_, N1,d_eff, LDA,sigma1, U_1, LDU, VT, LDVT,WORK, LWORK, RWORK, INFO )
    !*     Check for convergence.
    
    IF( INFO.GT.0 ) THEN
        WRITE(*,*)'The algorithm computing SVD failed to converge.'
        STOP
    END IF
    deallocate(WORK,RWORK)
    do i=1,M1_
        write(730,*) i,sigma1(i)
    enddo

    return
   
    
  end subroutine s_v_d