!--------------------------------------------------------------------------------
	module grid_data
	implicit none
	save
	integer::scheme_choice,n
	double precision::u_vel
	type scaler_node_data
		integer::bc
		double precision::x,y,temp
		double precision::fw,fe,dw,de
		double precision::aw,an,as,ae,ap,p_error,b
	end type
	type(scaler_node_data),dimension(100)::scaler_node
	double precision,dimension(100)::a,b,d,source,var
	end module grid_data
!------------------------------------------------------------------------------------
	program two_d_cd
	use grid_data
	implicit none
	integer::i,j
	n = 100
	scheme_choice = 1
	u_vel = 0.0001d0
	call compute_flux_and_diffusion_terms
	select case(scheme_choice)
		CASE(1) 
			call first_order_upwind_sch_coeff
		CASE(2)
			call central_sch_coeff
		CASE(3)
			call hybrid_differencing_scheme
		CASE default
	end select
	!do i=2,(n-1)
	 !scaler_node(i)%temp = 0.0d0
	!end do
	!scaler_node(1)%temp = 100.0d0
	!scaler_node(n)%temp = 0.0d0
	call tdma_initialization
	call tdma(n)
	do i=1,n
		scaler_node(i)%temp = var(i)
	end do
	call write_data	
	end program two_d_cd
!---------------------------------------------------------------------------------------
	subroutine compute_flux_and_diffusion_terms
	use grid_data
	implicit none
	integer::i,j
	do i=2,(n-1)
		scaler_node(i)%fe = 1000*u_vel        !Considering....density = 1000, u = 0.01m/s
		scaler_node(i)%fw = 1000*u_vel  
	end do
	do i=2,(n-1)
		scaler_node(i)%de = 111.1d0
		scaler_node(i)%dw = 111.1d0
	end do
	end subroutine compute_flux_and_diffusion_terms
!-------------------------------------------------------------------------------------------
	subroutine central_sch_coeff
	use grid_data
	implicit none
	integer::i
	do i=2,(n-1)
		scaler_node(i)%ae = scaler_node(i)%de - scaler_node(i)%fe/2.0d0  !max(-scaler_node(i)%fe,0.00)
		scaler_node(i)%aw = scaler_node(i)%dw + scaler_node(i)%fw/2.0d0  !max(scaler_node(i)%fw,0.00) 
		scaler_node(i)%ap = scaler_node(i)%ae + scaler_node(i)%aw + (scaler_node(i)%fe - scaler_node(i)%fw) 
	end do
	end subroutine central_sch_coeff
!-----------------------------------------------------------------------------------------
    subroutine first_order_upwind_sch_coeff
    use grid_data
    implicit none
    
	integer::i,j
         
	do i=2,(n-1)
		scaler_node(i)%ae = scaler_node(i)%de + max(-scaler_node(i)%fe,0.0d0)
		scaler_node(i)%aw = scaler_node(i)%dw + max(scaler_node(i)%fw,0.0d0) 
		scaler_node(i)%ap = scaler_node(i)%ae + scaler_node(i)%aw + (scaler_node(i)%fe - scaler_node(i)%fw) 
	end do
    end subroutine first_order_upwind_sch_coeff
!-----------------------------------------------------------------------------------------------------
	subroutine hybrid_differencing_scheme
	use grid_data
 	implicit none
	integer::i,j
	do i=2,(n-1)
		scaler_node(i)%ae = max(-scaler_node(i)%fe,(scaler_node(i)%de - (scaler_node(i)%fe)/2),0.0d0)
		scaler_node(i)%aw = max(scaler_node(i)%fw,(scaler_node(i)%dw + (scaler_node(i)%fw)/2),0.0d0) 
		scaler_node(i)%ap = scaler_node(i)%ae + scaler_node(i)%aw +(scaler_node(i)%fe - scaler_node(i)%fw) 
	end do
	end subroutine hybrid_differencing_scheme
!-----------------------------------------------------------------------------------------------------
	subroutine tdma_initialization
	use grid_data
	implicit none
	integer::i,j
	do i=2,(n-1)
		a(i) = scaler_node(i)%ae
		b(i) = scaler_node(i)%aw
		d(i) = scaler_node(i)%ap
		source(i) = 0.0d0
	end do
	d(1) = (3*111.1d0) + (1000*u_vel)/2.0d0	
	d(n) = (3*111.1d0) - (1000*u_vel)/2.0d0				
	b(1) = 0.0d0				
	a(1) = 111.1d0 - (1000*u_vel)/2.0d0
	b(n) = 111.1d0 + (1000*u_vel)/2.0d0				
	a(n) = 0.0d0			
	source(1) = (2*111.1d0 + 1000*u_vel)*(100.0d0)
	source(n) = 0.0d0
	end subroutine tdma_initialization
!-----------------------------------------------------------------------------------------------------
	subroutine tdma(nn)
	use grid_data
	implicit none
	integer::i,j
	integer,intent(in)::nn
	double precision::aa(nn),cc(nn)
	aa(1) = a(1)/(d(1) - 0.0d0)
	aa(nn) = 0.0d0							!- because a(nn) == 0
	do i=2,(nn-1)
		aa(i) = a(i)/(d(i) - (b(i)*aa(i-1)))
	end do
	cc(1) = (0.0d0 + source(1))/(d(1) - 0.0d0)
	do i=2,nn
		cc(i) = (b(i)*cc(i-1) + source(i))/(d(i) - (b(i)*aa(i-1)))
	end do
	var(nn) = 0.0d0 + cc(nn)
	do i=(nn-1),1,-1 
		var(i) = (aa(i)*var(i+1)) + cc(i)
	end do
	end subroutine tdma
!---------------------------------------------------------------------------------------------
	subroutine write_data
	use grid_data
	implicit none
	integer::i,j,done
	open(unit=15,file='temp.out',status='new',iostat=done)
        if(done/=0) then
			print*,"first remove the 'temp.out'file"
			stop
		end if
	do i=1,n
		write(15,*) scaler_node(i)%temp
	end do
	write(15,*) ""
	close(15)
	end subroutine write_data
!-------------------------------------------------------------------------------------------	
!-------------------------------------------------------------------------------------------------
