!------------------------------------------------------------------------------
	module grid_data
	implicit none
	save
	integer::scheme_choice,n
	double precision::u_vel,v_vel,fv,dv,fh,dh,ap0,top_bc,bottom_bc,right_bc,left_bc,heat_flux
	type scaler_node_data
		integer::bc
		double precision::x,y,temp,residual,temp_difference,temp_prev
		double precision::fw,fe,fn,fs,dw,de,dn,ds
		double precision::aw,an,as,ae,ap,p_error,b,su,sp
	end type
	type(scaler_node_data),dimension(100,100)::scaler_node
	double precision,dimension(100)::a,b,d,source,var
	double precision,dimension(100)::u  !-store the precalculated velocity field
	double precision,dimension(100)::v
	end module grid_data
!------------------------------------------------------------------------------
	program twod_fvm_cd
	use grid_data
	implicit none
	integer::i,j,xi,yi,ii,jj,conv_stat_outer,k,done
	n = 100
	scheme_choice = 2
	heat_flux = 10
	right_bc = 10 !- temp
	left_bc = 10 !- temp
	!del_y = 
	!del_x = 
	ap0 = 100.0d0
 
	call calculate_u_vel
	call calculate_v_vel
	call temperature_initialization
	call compute_flux_and_diffusion_terms
	call central_sch_coeff
	

	open(unit=15,file='temp.out',status='new',iostat=done)
        if(done/=0) then
			print*,"first remove the 'temp.out'file"
			stop
		end if

	do k=1,500000

		!print*, "Calculations start at time step", k
	  	call solve_temperature

		do i=1,n
		 do j=1,n
			scaler_node(i,j)%temp_difference = dabs(scaler_node(i,j)%temp - scaler_node(i,j)%temp_prev)
			scaler_node(i,j)%temp_prev = scaler_node(i,j)%temp
		 end do
		end do


		if((k - ((k/5000)*5000)) == 0) then
			call write_data
		end if


		conv_stat_outer = 0
		do i=1,n
		 do j=1,n
			if(scaler_node(i,j)%temp_difference>0.001) then
				conv_stat_outer = conv_stat_outer + 1
			end if
		 end do
		end do
		if(conv_stat_outer == 0) then
		  	 print*, "Solution coverged at time step",k
			 !call compute_average_values
			 call write_data
		  	 EXIT
			end if
		print*, "Calculations end at time step", k
	end do

	close(15)
		
	end program twod_fvm_cd

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
	subroutine calculate_u_vel
	use grid_data
	implicit none
	integer::i,j,xx,yy
	double precision::delta_y,h
	h = 0.01; !- metre
	delta_y = 0.0001; !- metre
	!- max velocity = 1.5623e-02 m/s
	!- velocity at the last two nodes = 3.1094e-04 m/s
	do j=1,n
		u(j) = -(0.01/(2*0.0008))*(((j-1)*delta_y + 0.00005 )*((j-1)*delta_y + 0.00005 ) - h*((j-1)*delta_y + 0.00005))
	end do
	
	end subroutine calculate_u_vel
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
	subroutine calculate_v_vel
	use grid_data
	implicit none
	integer::i,j,xx,yy
	do j=1,n
		v(j) = 0.0d0
	end do
	
	end subroutine calculate_v_vel
!--------------------------------------------------------------------------
	subroutine compute_flux_and_diffusion_terms
	use grid_data
	implicit none
	integer::i,j,xx,yy

	dh = 1.433!*del_y
	dv = 1.433!*del_x
	
	do xx = 1,n
		do yy = 1,n
			scaler_node(xx,yy)%fe = 1000*u(yy)
			scaler_node(xx,yy)%fw = 1000*u(yy)
			scaler_node(xx,yy)%fn = 1000*v(yy)
			scaler_node(xx,yy)%fs = 1000*v(yy)
		end do
	end do

	do xx = 1,n
		do yy = 1,n
			scaler_node(xx,yy)%de = dh !-can changed depending on delta_x and delta_y
			scaler_node(xx,yy)%dw = dh !-can changed depending on delta_x and delta_y
			scaler_node(xx,yy)%dn = dv !-can changed depending on delta_x and delta_y
			scaler_node(xx,yy)%ds = dv !-can changed depending on delta_x and delta_y
		end do
	end do


	end subroutine compute_flux_and_diffusion_terms

!-----------------------------------------------------------------------------
!- 2D coefficient settings ---------------------------------------------------

	subroutine central_sch_coeff
	use grid_data
	implicit none
	integer::i,j,xx,yy
	
	!- set the inner nodes first--------------------------------------------------
	do xx = 2,(n-1)
		do yy = 2,(n-1)
			scaler_node(xx,yy)%ae = dh - 1000*u(yy)/2.0d0
			scaler_node(xx,yy)%aw = dh + 1000*u(yy)/2.0d0
			scaler_node(xx,yy)%an = dv - 1000*v(yy)/2.0d0
			scaler_node(xx,yy)%as = dv + 1000*v(yy)/2.0d0
			scaler_node(xx,yy)%su = 0.0d0	!- inner nodes
			scaler_node(xx,yy)%sp = 0.0d0	!- inner nodes
			!- now sum up all the terms to get the aP along with flux terms and -sp
			scaler_node(xx,yy)%ap = scaler_node(xx,yy)%ae + scaler_node(xx,yy)%aw + scaler_node(xx,yy)%an &
						&+ scaler_node(xx,yy)%as + ap0 - scaler_node(xx,yy)%sp 
		end do
	end do

	!-//--------------------------------------------------------------------------
	!- set TOP and BOTTOM sides excluding corners---------------------------------------
	do xx = 2,(n-1)
		!- top face---------CV 7 -----excluding corner----------------------------
		!- const heat flux = 100 KW/m.m
		scaler_node(xx,n)%ae = dh - 1000*u(n)/2.0d0
		scaler_node(xx,n)%aw = dh + 1000*u(n)/2.0d0
		scaler_node(xx,n)%an = 0.0d0 				!-dv - 1000*v(n)
		scaler_node(xx,n)%as = dv + 1000*v(n)/2.0d0
		scaler_node(xx,n)%su = heat_flux*1000*0.01 		!-(2*dv - fv)*top_bc ! here dn is double fn is same
		scaler_node(xx,n)%sp = 1000*v(n)  		! here dn is double and fn is same
		scaler_node(xx,n)%ap = scaler_node(xx,n)%ae + scaler_node(xx,n)%aw + scaler_node(xx,n)%an + &
					& scaler_node(xx,n)%as + ap0 - scaler_node(xx,n)%sp 
		
		!- bottom face---------CV 5--------------excluding corner---------------------- 
		!- const heat flux = 100 KW/m.m
		scaler_node(xx,1)%ae = dh - 1000*u(1)/2.0d0
		scaler_node(xx,1)%aw = dh + 1000*u(1)/2.0d0
		scaler_node(xx,1)%an = dv - 1000*v(1)/2.0d0
		scaler_node(xx,1)%as = 0.0d0
		scaler_node(xx,1)%su = heat_flux*1000*0.01 ! here ds is double fs is same
		scaler_node(xx,1)%sp = -1000*v(1)         ! here ds is double and fs is same
		scaler_node(xx,1)%ap = scaler_node(xx,1)%ae + scaler_node(xx,1)%aw + scaler_node(xx,1)%an +&
					& scaler_node(xx,1)%as + ap0 - scaler_node(xx,1)%sp 
		
	end do
	
	!-//--------------------------------------------------------------------------
	!- set LEFT and RIGHT sides excluding corners---------------------------------------------------

	do yy = 2,(n-1)
		!- Left face---------CV 8 -----excluding corner----------------------------
		scaler_node(1,yy)%ae = dh - 1000*u(yy)/2.0d0
		scaler_node(1,yy)%aw = 0.0d0 !ok
		scaler_node(1,yy)%an = dv - 1000*v(yy)/2.0d0
		scaler_node(1,yy)%as = dv + 1000*v(yy)/2.0d0
		scaler_node(1,yy)%su = (1000*u(yy) + 2*dh)*left_bc ! here dw is double fw is same
		scaler_node(1,yy)%sp = -(1000*u(yy) + 2*dh) !- here dw is double and fw is same
		scaler_node(1,yy)%ap = scaler_node(1,yy)%ae + scaler_node(1,yy)%aw + scaler_node(1,yy)%an +&
					& scaler_node(1,yy)%as + ap0 - scaler_node(1,yy)%sp
		
		!- Right face---------CV 6--------------excluding corner---------------------- 
		scaler_node(n,yy)%ae = 0.0d0 !ok
		scaler_node(n,yy)%aw = dh + u(yy)/2.0d0
		scaler_node(n,yy)%an = dv - v(yy)/2.0d0
		scaler_node(n,yy)%as = dv + v(yy)/2.0d0
		scaler_node(n,yy)%su = -(1000*u(yy) - 2*dh)*right_bc ! here de is double fe is same
		scaler_node(n,yy)%sp = -(1000*u(yy) + 2*dh) ! here de is double and fe is same
		scaler_node(n,yy)%ap = scaler_node(n,yy)%ae + scaler_node(n,yy)%aw + scaler_node(n,yy)%an +&
					& scaler_node(n,yy)%as + ap0 - scaler_node(n,yy)%sp
		
	end do
	!-//--------------------------------------------------------------------------
	!- corners co-efficients -----------------------------------------------------

		
	!-------------CV 4 -----top left corner--------------------------------------

		scaler_node(1,n)%ae = dh - 1000*u(n)/2.0d0
		scaler_node(1,n)%aw = 0.0d0
		scaler_node(1,n)%an = 0.0d0
		scaler_node(1,n)%as = dv + 1000*v(n)/2.0d0
		scaler_node(1,n)%su = (1000*u(yy) + 2*dh)*left_bc + heat_flux*1000*0.01 !- here dn is double fn is same
		scaler_node(1,n)%sp = -(1000*u(yy) + 2*dh) + 1000*v(n)	  ! here dn is double and fn is same
		scaler_node(1,n)%ap = scaler_node(1,n)%ae + scaler_node(1,n)%aw + scaler_node(1,n)%an + &
					&scaler_node(1,n)%as + ap0 - scaler_node(1,n)%sp
		

	!-------------CV 3 -----top right corner--------------------------------------


		scaler_node(n,n)%ae = 0.0d0
		scaler_node(n,n)%aw = dh + 1000*u(n)/2.0d0
		scaler_node(n,n)%an = 0.0d0
		scaler_node(n,n)%as = dv + 1000*v(n)/2.0d0
		scaler_node(n,n)%su = heat_flux*1000*0.01 -(1000*u(yy) - 2*dh)*right_bc ! here dn is double fn is same
		scaler_node(n,n)%sp = -(1000*u(yy) + 2*dh) + 1000*v(n) 	  ! here dn is double and fn is same
		scaler_node(n,n)%ap = scaler_node(n,n)%ae + scaler_node(n,n)%aw + scaler_node(n,n)%an + &
					&scaler_node(n,n)%as + ap0 - scaler_node(n,n)%sp
		


	!---------------CV 1 ------bottom left corner----------------------------------------

		scaler_node(1,1)%ae = dh - 1000*u(1)/2.0d0
		scaler_node(1,1)%aw = 0.0d0
		scaler_node(1,1)%an = dv - 1000*v(1)/2.0d0
		scaler_node(1,1)%as = 0.0d0
		scaler_node(1,1)%su = heat_flux*1000*0.01 + (1000*u(yy) + 2*dh)*left_bc
		scaler_node(1,1)%sp = -1000*v(1) -(1000*u(yy) + 2*dh)
		scaler_node(1,1)%ap = scaler_node(1,1)%ae + scaler_node(1,1)%aw + scaler_node(1,1)%an +&
					& scaler_node(1,1)%as + ap0 - scaler_node(1,1)%sp

	!-----------------CV 2------------ bottom right corner------------------------------------

		scaler_node(n,1)%ae = 0.0d0
		scaler_node(n,1)%aw = dh + 1000*u(1)/2.0d0
		scaler_node(n,1)%an = dv - 1000*v(1)/2.0d0
		scaler_node(n,1)%as = 0.0d0
		scaler_node(n,1)%su = heat_flux*1000*0.01 - (1000*u(yy) - 2*dh)*right_bc
		scaler_node(n,1)%sp = -(1000*u(yy) + 2*dh) - 1000*v(1)  
		scaler_node(n,1)%ap = scaler_node(n,1)%ae + scaler_node(n,1)%aw + scaler_node(n,1)%an +&
					& scaler_node(n,1)%as + ap0 - scaler_node(n,1)%sp

!--------------------------------------------------------------------------------------------------------------------


	end subroutine central_sch_coeff


!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------

	subroutine solve_temperature
	use grid_data
	implicit none

	integer::i,j,xi,yi,ii,jj,conv_stat


	
	!-//sweep loop starts----------------------------------------

	do i=1,10000

		!print*, "                          Iteration", i, "starts"


			
	 		do xi = 1,n 	!- sweep along x-direction

		 	    call tdma_initialization(xi)
		 	    call tdma(n)	
		
			    do yi = 1,n
				scaler_node(xi,yi)%residual = dabs(var(yi) - scaler_node(xi,yi)%temp)
				scaler_node(xi,yi)%temp = var(yi)
		 	    end do

	  		end do

			!-//sweep loop ends----------------------------------------

			conv_stat = 0
			do ii=1,n
		  	 do jj=1,n
				if(scaler_node(ii,jj)%residual>0.0001) then
				conv_stat = conv_stat + 1
				end if
		   	 end do
			end do
			if(conv_stat == 0) then
		  	 !print*, "Solution coverged at interation number",i
		  	 EXIT
			end if

			!print*, "                          Iteration", i, "starts"
		end do




	end subroutine solve_temperature


!-----------------------------------------------------------------------------------------------------

	SUBROUTINE temperature_initialization
	use grid_data
	implicit none

	integer::i,j

	do i=1,n
	 do j=1,n
	  scaler_node(i,j)%temp = 10.0d0
	 end do
	end do

	do i=1,n
	 do j=1,n
		scaler_node(i,j)%temp_prev = 10.0d0
	 end do
	end do


	end subroutine temperature_initialization

!-------------------------------------------------------------------------------------------------------------


	!-----------------------------------------------------------------------------------------------------
	!- tdma initialization for each and every row/ or column to be processed before sending to tdma(xi)
	!- definition od actual source terms and tridiagonal structure creation with b(i),d(i),and a(i)
	!- where -a(i) = aE(i) d(i) = aP(i) and -b(i) = aW(i)--------------------------------------------------

	subroutine tdma_initialization(col)
	use grid_data
	implicit none

	integer::i,j
	integer,intent(in)::col


	do i=1,n
		b(i) = scaler_node(col,i)%as
	end do


	do i=1,n
		a(i) = scaler_node(col,i)%an
	end do

	do i=1,n
		d(i) = scaler_node(col,i)%ap
	end do

	do i=1,n
		if(col==1) then
		 source(i) = scaler_node(col,i)%su + scaler_node(col,i)%ae*scaler_node((col+1),i)%temp + ap0*scaler_node(col,i)%temp_prev	
		else if(col==n) then
		 source(i) = scaler_node(col,i)%su + scaler_node(col,i)%aw*scaler_node((col-1),i)%temp + ap0*scaler_node(col,i)%temp_prev
		else if((col/=1).and.(col/=n)) then
		 source(i) = scaler_node(col,i)%su + scaler_node(col,i)%aw*scaler_node((col-1),i)%temp &
				&+ scaler_node(col,i)%ae*scaler_node((col+1),i)%temp + ap0*scaler_node(col,i)%temp_prev
		end if
	end do



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


!-------------------------------------------------------------------------------------------------------



!--------------------------------------------------------------------------------------------------------

	subroutine write_data
	use grid_data
	implicit none
	integer::i,j,done
	
	do i=1,n
		do j=1,n
			write(15,*) scaler_node(i,j)%temp
		end do
	 	!--write(15,*) "...........line", i, "finished.........."	
	end do
	!write(15,*) ""
	


	end subroutine write_data
	!-------------------------------------------------------------------------------------------	
	!-------------------------------------------------------------------------------------------------


