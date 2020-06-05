SUBROUTINE despekl(maxgate,maxazim,                                     &
                   ngate,nazim,varchek,medfilt,rdata,tmp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Despeckle radar data.
!  Can be used for reflectivity or velocity.
!
!
!  AUTHOR: Keith Brewster
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    maxgate  Maximum number of gates in a radial
!    maxazim  Maximum number of radials in a tilt
!    ngate    Number of gates in radial
!    nazim    Number of radials
!    varchek  Threshold to determine good data vs. flagged
!    medfilt  Median filter option switch
!    rdata    Doppler radial velocity
!
!  OUTPUT:
!    rdata
!
!  WORK ARRAYS:
!
!    tmp
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: maxgate
  INTEGER, INTENT(IN) :: maxazim
  INTEGER, INTENT(IN) :: ngate
  INTEGER, INTENT(IN) :: nazim

  REAL,    INTENT(IN) :: varchek
  INTEGER, INTENT(IN) :: medfilt

  REAL, INTENT(INOUT) :: rdata(maxgate,maxazim)
  REAL, INTENT(OUT)   :: tmp(maxgate,maxazim)

!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  INTEGER :: NN,m,n,loc,is
  INTEGER :: kntgd,kntdsp
  INTEGER :: nprint
!
  REAL :: sortdata(9)
  REAL :: swp
!
  REAL :: sum,pctdsp
  REAL, PARAMETER :: dspmiss = -991.
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
! Create temporary array where tmp=1 if non-missing tmp=0 if missing
!
!-----------------------------------------------------------------------
!
!
  tmp=0.
  DO j=1,nazim
    DO i=1,ngate
      IF ( rdata(i,j) > varchek ) tmp(i,j)=1.
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
! If datum has fewer than 3 non-missing neighbors in a 3x3
! square template set it to missing
!
!-----------------------------------------------------------------------
!
  kntgd=0
  kntdsp=0
  DO j=2,nazim-1
    DO i=2,ngate-1
      IF(rdata(i,j) > varchek ) THEN
        kntgd=kntgd+1
        sum=tmp(i-1,j+1)+tmp(i,j+1)+tmp(i+1,j+1)          &
           +tmp(i-1,j  )+           tmp(i+1,j)            &
           +tmp(i-1,j-1)+tmp(i,j-1)+tmp(i+1,j-1)
        IF( sum < 3. ) THEN
          kntdsp=kntdsp+1
          rdata(i,j) = dspmiss
        END IF
      END IF
    END DO
  END DO
  IF(kntgd > 0 ) THEN
    pctdsp=100.*(float(kntdsp)/float(kntgd))
  ELSE
    pctdsp=0.
  END IF

  write(6,'(a,i8,a,i8,a,f6.1,a)') &
    ' Despeckled ',kntdsp,' of ',kntgd,' data =',pctdsp,' percent.'

!----------------------------------------------------------------------------
!
!  Zhao Kun added a median filter to smooth the reflectivity and velocity
!  Recoded by Keith Brewster to make it optional and streamline a bit.
!
!----------------------------------------------------------------------------

  IF(medfilt > 0) THEN
!   nprint=0
    DO j=1,nazim
      DO i=1,ngate
        tmp(i,j)=rdata(i,j)
      END DO
    END DO

    DO j=2,nazim-1
      DO i=2,ngate-1
        NN=0
        DO m=-1,1
          DO n=-1,1
            swp=tmp(i+m,j+n)
            IF(swp > varchek) THEN
              NN=NN+1
              IF(NN == 1) THEN ! first
                sortdata(NN)=swp
              ELSE
                DO loc=1,NN-1
                  IF(swp < sortdata(loc)) THEN
                    DO is=NN,loc+1,-1
                      sortdata(is)=sortdata(is-1)
                    END DO
                    EXIT
                  END IF
                END DO
                sortdata(loc)=swp
              END IF
            END IF
          END DO
        END DO

        IF (NN > 6) THEN
!         IF(nprint < 100) THEN
!           print *, ' NN: ',NN
!           print *, ' sortdata: ',(sortdata(m),m=1,NN)
!           print*,'median',sortdata((NN+1)/2),rdata(i,j),tmp(i,j)
!           nprint=nprint+1
!         END IF
          rdata(i,j)=sortdata((NN+1)/2)

        END IF

      END DO
    END DO
!
!  First radial
!
    j=1
    DO i=2,ngate-1
      NN=0
      DO m=-1,1
        DO n=0,2
          swp=tmp(i+m,j+n)
          IF(swp > varchek) THEN
            NN=NN+1
            IF(NN == 1) THEN ! first
              sortdata(NN)=swp
            ELSE
              DO loc=1,NN-1
                IF(swp < sortdata(loc)) THEN
                  DO is=NN,loc+1,-1
                    sortdata(is)=sortdata(is-1)
                  END DO
                  EXIT
                END IF
              END DO
              sortdata(loc)=swp
            END IF
          END IF
        END DO
      END DO

      IF (NN > 6) THEN

!       DO m=1,NN
!         print*,'NN',m,sortdata(m)
!       ENDDO
!       print*,'mean',sortdata((NN+1)/2),rdata(i,j),tmp(i,j)
        rdata(i,j)=sortdata((NN+1)/2)
!       stop

      END IF

    END DO
!
!  Last radial
!
    j=nazim
    DO i=2,ngate-1
      NN=0
      DO m=-1,1
        DO n=-2,0
          swp=tmp(i+m,j+n)
          IF(swp > varchek) THEN
            NN=NN+1
            IF(NN == 1) THEN ! first
              sortdata(NN)=swp
            ELSE
              DO loc=1,NN-1
                IF(swp < sortdata(loc)) THEN
                  DO is=NN,loc+1,-1
                    sortdata(is)=sortdata(is-1)
                  END DO
                  EXIT
                END IF
              END DO
              sortdata(loc)=swp
            END IF
          END IF
        END DO
      END DO

      IF (NN > 6) THEN

!       DO m=1,NN
!         print*,'NN',m,sortdata(m)
!       ENDDO
!       print*,'mean',sortdata((NN+1)/2),rdata(i,j),tmp(i,j)
        rdata(i,j)=sortdata((NN+1)/2)
!       stop

      END IF

    END DO
  END IF ! median filter

  RETURN
!
END SUBROUTINE despekl
