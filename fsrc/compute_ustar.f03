!***********************************************************************
subroutine compute_ustar(distanza,reynol,umedia,u_t,rough,rougheight)
    !***********************************************************************
    !
    implicit none
    !
    !-----------------------------------------------------------------------
    !     variables declaration
    real transition,rougheight,kdynamic,coef_rough
    real errore_star,reynol,z0
    real umedia,u_t,distanza
    real argomentolog,f,fprime,tempv
    integer contatore_star,i,rough
    !
    !-----------------------------------------------------------------------
    ! parameter  log profile

    transition = 0
    contatore_star = 0
    kdynamic = 2.44   ! 1/k with k von karman constant
    coef_rough = 5.1
    errore_star = 1000.
    !-----------------------------------------------------------------------
    !
    !     CASE 1: smooth surface

    do i=1,1-rough
        !     Newton - Rapshon iterative procedure
        do while(errore_star>1.d-4 .and. contatore_star<100)
            tempv=u_t
            u_t = abs(u_t)
         	 
            argomentolog = 1.+ u_t*rougheight*reynol
	 
            coef_rough = coef_rough  &
                - real(transition)*kdynamic*log(argomentolog)

            f=u_t*(kdynamic*log(abs(distanza*u_t*reynol))+ coef_rough) - umedia
     
            !           fprime is f derivative
            fprime = kdynamic*log(abs(distanza*u_t*reynol))+kdynamic+coef_rough &
            -real(transition)*u_t*kdynamic*rougheight/(1./reynol + rougheight*u_t)
	 
            u_t = u_t - f/fprime

            if(u_t>1.E-8)then
                errore_star=abs(u_t-tempv)
            !v	     errore_star = abs(f/sqrt(u_t))
            end if
	 
            contatore_star = contatore_star +1
            coef_rough = 5.1

        end do   !end loop do while
    end do !smooth surface

    !-----------------------------------------------------------------------
    !
    !     CASE 2: rough surface

    do i=1,rough
        u_t=0.41*umedia/log(distanza/rougheight)
    end do
      
      
    return
end
