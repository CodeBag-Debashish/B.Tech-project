!--------------------------------------------------------------------------------

	module grid_data
	implicit none
	save


	integer::scheme_choice,n, temp_left, temp_right
	double precision::ap0

	type scaler_node_data
         integer::bc
         double precision::x,y,temp,temp_residual
	 double precision::fw,fe,dw,de
         double precision::aw,an,as,ae,ap,b,su,sp
    	end type

	type(scaler_node_data),dimension(100)::scaler_node
	double precision,dimension(100)::a,b,d,source,var


	end module grid_data


!------------------------------------------------------------------------------------



	program one_d_cd_unsteady
	use grid_data
	implicit none


	integer::i,j,k,conv_stat,done

	n = 100
	scheme_choice = 2
	ap0 = 1000.0d0

	temp_left = 100.0d0
	temp_right = 0.0d0

	do i=1,n
	  scaler_node(i)%temp = 0.0d0
	end do

	call compute_flux_and_diffusion_terms
	call central_sch_coeff


	open(unit=15,file='temp.out',status='new',iostat=done)
        if(done/=0) then
      	 print*,"first remove the 'temp.out'file"
      	stop
     	end if

	do k=1,500000

	print*, "Calculations start for time step", k

	

	 	call tdma_initialization
		call tdma(n)

		do i=1,n
		 scaler_node(i)%temp_residual = dabs(scaler_node(i)%temp - var(i))
	 	 scaler_node(i)%temp = var(i)
		end do

		if((k - ((k/100)*100)) == 0) then
			call write_data
		end if


		conv_stat = 0
		do i=1,n
		 if(scaler_node(i)%temp_residual>0.00001) then
			conv_stat = conv_stat + 1
		 end if
		end do


		if(conv_stat==0) then
			print*, "Solution converged"
			call write_data
			STOP
		end if
	print*, "Calculations end for time step", k


	end do

	close(15)
	

	end program one_d_cd_unsteady



!---------------------------------------------------------------------------------------



	subroutine compute_flux_and_diffusion_terms
	use grid_data
	implicit none

	integer::i,j


	do i=1,n
		scaler_node(i)%fe = 1.0d0!1000*0.01       !Considering....density = 1000, u = 0.0001m/s
		scaler_node(i)%fw = 1.0d0!1000*0.01   
	end do


	do i=1,n
		scaler_node(i)%de = 111.1d0	!Considering k=0.6, cp=4186J/kgK, del_x=0.001
		scaler_node(i)%dw = 111.1d0
	end do
	
	scaler_node(n)%de = 2*111.1d0
	scaler_node(1)%dw = 2*111.1d0


	end subroutine compute_flux_and_diffusion_terms


!-------------------------------------------------------------------------------------------

    subroutine central_sch_coeff
    use grid_data
    implicit none
    
    integer::i
         

     do i=1,n
       scaler_node(i)%ae = scaler_node(i)%de - scaler_node(i)%fe/2.0d0  !max(-scaler_node(i)%fe,0.00)
       scaler_node(i)%aw = scaler_node(i)%dw + scaler_node(i)%fw/2.0d0  !max(scaler_node(i)%fw,0.00) 
	scaler_node(i)%su = 0.0d0
	scaler_node(i)%sp = 0.0d0 

       scaler_node(i)%ap = scaler_node(i)%ae + scaler_node(i)%aw + (scaler_node(i)%fe - scaler_node(i)%fw) &
				&+ ap0 - scaler_node(i)%sp
     end do


	scaler_node(1)%aw = 0.0d0
	scaler_node(1)%sp = -(scaler_node(1)%dw + scaler_node(1)%fw)
	scaler_node(1)%su = (scaler_node(1)%dw + scaler_node(1)%fw)*temp_left
	scaler_node(1)%ap = scaler_node(1)%ae + scaler_node(1)%aw + (scaler_node(1)%fe - scaler_node(1)%fw) &
				&+ ap0 - scaler_node(1)%sp
	scaler_node(n)%ae = 0.0d0
	scaler_node(n)%sp = scaler_node(n)%fe - scaler_node(n)%de
	scaler_node(n)%su = (scaler_node(n)%de - scaler_node(n)%fe)*temp_right
	scaler_node(n)%ap = scaler_node(n)%ae + scaler_node(n)%aw + (scaler_node(n)%fe - scaler_node(n)%fw) &
				&+ ap0 - scaler_node(n)%sp

 

   end subroutine central_sch_coeff

!-----------------------------------------------------------------------------------------



	subroutine tdma_initialization
	use grid_data
	implicit none

	integer::i,j
	

	do i=1,n
	 a(i) = scaler_node(i)%ae
	 b(i) = scaler_node(i)%aw
	 d(i) = scaler_node(i)%ap
	 source(i) = scaler_node(i)%su + ap0*scaler_node(i)%temp
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


!---------------------------------------------------------------------------------------------


	subroutine write_data
	use grid_data
	implicit none
	
	integer::i,j,done
	

	

	do i=1,n
	   write(15,*) scaler_node(i)%temp
	end do
	
	!write(15,*) ""





	end subroutine write_data

	

!-------------------------------------------------------------------------------------------	


!-------------------------------------------------------------------------------------------------
