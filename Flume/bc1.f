
c
c
c     =================================================================
      subroutine bc1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt,mthbc)
c     =================================================================
c
c     # Standard boundary condition choices for claw2
c
c     # Modified for 1d GeoClaw to extend topo aux(1,:) to ghost cells.
c
c     # At each boundary  k = 1 (left),  2 (right):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd component of q.
c     ------------------------------------------------
c
c     # Extend the data from the computational region
c     #      i = 1, 2, ..., mx2
c     # to the virtual cells outside the region, with
c     #      i = 1-ibc  and   i = mx+ibc   for ibc=1,...,mbc
c
      use geoclaw_module, only: grav
      use myfunctions, only: racine3, root3
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)

      dimension mthbc(2)
      dimension rootH(3)

      pi = 4.d0*atan(1.d0)
      h0 = 2.d0
      q0 = 10
      p = 0.5 ! pelle
      
      ampl_eta = 0.03
      ampl_u = ampl_eta * sqrt(grav/h0)
      !write(6,*) '+++ ampl_eta, ampl_u: ',ampl_eta, ampl_u
      !ampl = 0.047  #0.15d0
      tperiod = 20.d0

c
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      do ibc=1,mbc
         aux(1,1-ibc)=aux(1,1)
      enddo

      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
c     # wave maker
c      if (t < tperiod) then
c          s = ampl_u * sin(2.d0*pi*t/tperiod)
c      else
c          s = 0.d0
c      endif
      hur = q(2,1)
      hr = q(1,1)
      ur = hur/hr
      r = ur-2.0d0*sqrt(grav*hr)
     
      do ibc=1,mbc
         aux(1,1-ibc) = aux(1,ibc)
         if (t < 5) then
            q(1,1-ibc) = h0
            q(2,1-ibc) = q0
         else if (t < 10) then
            q(1,1-ibc) = h0*(10-t)/5
            q(2,1-ibc) = q0*((10-t)/5)**2
         else
            q(1,1-ibc) = 0
            q(2,1-ibc) = 0
         end if
      enddo

      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,1)
         end do
      end do
      go to 199

  120 continue
c     # periodic:
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,mx+1-ibc)
         end do
      end do
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,ibc)
         end do
c        # negate the normal velocity:
         q(2,1-ibc) = -q(2,ibc)
      end do
      go to 199

  199 continue

c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      do ibc=1,mbc
         aux(1,mx+ibc)=aux(1,mx)
      enddo

      go to (200,210,220,230) mthbc(2)+1
c
  200 continue
c     # user-specified boundary conditions go here in place of error output
c      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
c      stop
c      go to 299
      
c      do ibc=1,mbc
c         do m=2,meqn
c            q(m,mx+ibc) = q(m,mx)
c         end do
c         q(1,mx+ibc) = 1.d0
c      end do
c      go to 299
      hm = q(1,mx)
      qm = q(2,mx)
      
      if (hm > p) then
        Hb = p+1.5d0*(qm**2/grav)**(1.d0/3.d0)
        rootH = root3(-Hb,0.d0,qm**2/2/grav)
        hp = rootH(2)
        do ibc=1,mbc
             q(1,mx+ibc) = hp
             q(2,mx+ibc) = qm
        end do
      endif
      if (hm <= p) then
        do ibc=1,mbc
             q(1,mx+ibc) = hm
             q(2,mx+ibc) = 0.d0
        end do
      endif
      go to 299

  210 continue
c     # zero-order extrapolation:
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,mx)
         end do
      end do
      go to 299

  220 continue
c     # periodic:
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,ibc)
         end do
      end do
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,mx+1-ibc)
         end do
c        # negate the normal velocity:
         q(2,mx+ibc) = -q(2,mx+1-ibc)
      end do
      go to 299

  299 continue
c
      return
      end
