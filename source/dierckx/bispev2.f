      subroutine bispev2(tx,nx,ty,ny,c,kx,ky,x,y,mxy,z,ier)
c  Modified by D.Kelson (2003) to do arbitrary array of locations.
c
c  ..scalar arguments..
      integer nx,ny,kx,ky,mxy,lwrk,kwrk,ier
c  ..array arguments..
      real*8 tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mxy),y(mxy),z(mxy)
c  ..local scalars..
      integer i,iw,lwest
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 0
      call fpbisp2(tx,nx,ty,ny,c,kx,ky,x,y,mxy,z)
 100  return
      end
