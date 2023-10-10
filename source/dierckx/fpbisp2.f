      subroutine fpbisp2(tx,nx,ty,ny,c,kx,ky,x,y,mxy,z)
c       Attempt by D. Kelson (2003) to do the 2D b-spline evaluation
c       at array of arbitrary locations.
c  ..scalar arguments..
      integer nx,ny,kx,ky,mxy
c  ..array arguments..
      integer lc,lc1,lc2
      real*8 tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mxy),y(mxy),z(mxy)
c  ..local scalars..
      integer kx1,ky1,xl,xl1,xl2,yl,yl1,yl2,m,nkx1,nky1
      real*8 argx,argy,sp,txb,txe,tyb,tye
c  ..local arrays..
      real*8 hx(6),hy(6)
c  ..subroutine references..
c    fpbspl
c  ..
      kx1 = kx+1
      nkx1 = nx-kx1
      txb = tx(kx1)
      txe = tx(nkx1+1)

      ky1 = ky+1
      nky1 = ny-ky1
      tyb = ty(ky1)
      tye = ty(nky1+1)

      do 40 i=1,mxy
        l = kx1
        l1 = l+1
        argx = x(i)
        if(argx.lt.txb) argx = txb
        if(argx.gt.txe) argx = txe
c        do 401 l=kx1,nkx1
c          if (argx.ge.tx(l+1)) goto 20
c 401    continue
  10    if(argx.lt.tx(l1) .or. l.eq.nkx1) go to 20
        l = l1
        l1 = l+1
        go to 10
  20    call fpbspl(tx,nx,kx,argx,l,hx)
        xl = l-kx1

        l = ky1
        l1 = l+1
        argy = y(i)
        if(argy.lt.tyb) argy = tyb
        if(argy.gt.tye) argy = tye
c        do 501 l=ky1,nky1
c          if (argy.ge.ty(l+1)) goto 60
c 501    continue
  50    if(argy.lt.ty(l1) .or. l.eq.nky1) go to 60
        l = l1
        l1 = l+1
        go to 50
  60    call fpbspl(ty,ny,ky,argy,l,hy)
        yl = l-ky1

        l = xl*nky1
        l1 = l + yl
        sp = 0.
        do 110 i1=1,kx1
          l2 = l1
          do 100 j1=1,ky1
            l2 = l2 + 1
            sp = sp+c(l2)*hx(i1)*hy(j1)
 100      continue
          l1 = l1 + nky1
 110    continue
        z(i) = sp
  40  continue

      return
      end

