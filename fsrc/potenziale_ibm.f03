!***********************************************************************
subroutine potenziale_ibm(myid,nproc,tipo,deepl,deepr,kparasta, &
    kparaend,rightpe,leftpe,tagls,taglr,tagrs, &
    tagrr)
      
    !***********************************************************************
    ! pressure correction inside Immersed Boundary is setted equal to that
    ! outside
    !-----------------------------------------------------------------------
    use myarrays_velo3
    !
    use scala3
    use tipologia
    !
    use mpi

    implicit none

    !-----------------------------------------------------------------------
    !     array declaration
    integer i,j,k
    integer nproc,myid
    integer kparasta,kparaend
    integer req1,req2,req3,req4
    integer rightpe ,leftpe
    integer rightpem,leftpem
    integer ierr,istatus,status(MPI_STATUS_SIZE)
    integer tagls,taglr,tagrs,tagrr
    integer conto_ib

    integer deepl,deepr
    integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
      
    real fi_old
    !-----------------------------------------------------------------------
    !     solid cell  fi(solida)=fi(ib)

    do k=kparasta,kparaend
        do j=1,jy
            do i=1,jx
   
                if(tipo(i,j,k).eq.0)then
                    fi_old = fi(i,j,k)
                    fi(i,j,k)=0.
                    conto_ib=0
	
                    if(tipo(i+1,j  ,k  ).eq.1)then
                        fi(i  ,j  ,k  )=fi(i,j,k)+fi(i+1,j  ,k  )
                        conto_ib=conto_ib+1
                    end if
	  	  
                    if(tipo(i-1,j  ,k  ).eq.1)then
                        fi(i  ,j  ,k  )=fi(i,j,k)+fi(i-1,j  ,k  )
                        conto_ib=conto_ib+1
                    end if
	  	    
                    if(tipo(i  ,j+1,k  ).eq.1)then
                        fi(i  ,j  ,k  )=fi(i,j,k)+fi(i  ,j+1,k  )
                        conto_ib=conto_ib+1
                    end if
	  	    
                    if(tipo(i  ,j-1,k  ).eq.1)then
                        fi(i  ,j  ,k  )=fi(i,j,k)+fi(i  ,j-1,k  )
                        conto_ib=conto_ib+1
                    end if
	  	    
                    if(tipo(i  ,j  ,k+1).eq.1)then
                        fi(i  ,j  ,k  )=fi(i,j,k)+fi(i  ,j  ,k+1)
                        conto_ib=conto_ib+1
                    end if
	  	    
                    if(tipo(i  ,j  ,k-1).eq.1)then
                        fi(i  ,j  ,k  )=fi(i,j,k)+fi(i  ,j  ,k-1)
                        conto_ib=conto_ib+1
                    end if
	  
                    if(conto_ib.ne.0)then
                        fi(i,j,k)=fi(i,j,k)/conto_ib
                    else
                        fi(i,j,k)=fi_old
                    end if
	    
                end if
     
            end do
        end do
    end do

    !     send left
    if(myid.ne.0)then
        call MPI_SSEND(fi(0,0,kparasta),(jx+2)*(jy+2),MPI_REAL_SD, &
            leftpe,tagls,MPI_COMM_WORLD,ierr)
    !      call MPI_WAIT (req1,istatus,ierr)
    end if
    if(myid.ne.nproc-1)then
        call MPI_RECV(fi(0,0,kparaend+1),(jx+2)*(jy+2),MPI_REAL_SD, &
            rightpe,tagrr,MPI_COMM_WORLD,status,ierr)
    !      call MPI_WAIT (req2,istatus,ierr)
    end if

    !     send right
    if(myid.ne.nproc-1)then
        call MPI_SSEND(fi(0,0,kparaend),(jx+2)*(jy+2),MPI_REAL_SD, &
            rightpe,tagrs,MPI_COMM_WORLD,ierr)
    !      call MPI_WAIT (req3,istatus,ierr)
    end if
    if(myid.ne.0)then
        call MPI_RECV(fi(0,0,kparasta-1),(jx+2)*(jy+2),MPI_REAL_SD, &
            leftpe,taglr,MPI_COMM_WORLD,status,ierr)
    !      call MPI_WAIT (req4,istatus,ierr)
    end if
      
    return
end
