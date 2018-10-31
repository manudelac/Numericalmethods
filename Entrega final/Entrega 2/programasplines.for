c Manuel de la Cruz González 70909708H
      program programasplines

      real*8 x,f,xi,xf,width,y,z(0:500),s(0:500),j,a,b,c,d
      dimension a(0:500),d(0:500),c(0:500),b(0:500),f(0:500)
      dimension x(0:500)
      integer i,n1,p,k
      n1=21
c En este bucle leemos los datos del problema     
      open (55,file='splines.txt')
        do i=0,20
         read (55,*)x(i),f(i)
	    end do
      close (33)

c Pedimos los valores necesarios      
      write(*,*)'Introduzca el punto inicial del intervalo'
      read(*,*)xi
      write(*,*)'Introduzca el punto final del intervalo'
      read(*,*)xf
      write(*,*)'Introduce el numero de puntos'
      read(*,*)puntos
      
      width=(xf-xi)/(puntos)
      y=(xf-xi)/(n1-1)

c este do crea x equiespaciadas      
	  do i=0, n1-1
         x(i)=xi+y*i
      end do    

c con esta subrutina calcularemos los coeficientes     
      call calculosplines (n1,xi,xf,f,a,b,c,d,x)
c este vector son los puntos donde interpolamos     
	  do i=0, puntos-1
        z(i)=xi+width*i
      end do 

      do i=0, puntos-1
        
      open(30,file='interpolacion.txt')
c excepcion de que un punto sea extremo superior      
       if (z(i).ne.1.0d0) then
           j=(z(i)-xi)/y
           k=int(j)
        else
           k=n1-2
       end if
c en s estan las interpolaciones de z. Escribimos los puntos y la interpolacion en el fichero interpolacion       
      s(i)=a(k)*(z(i)-x(k))**3+b(k)*(z(i)-x(k))**2+c(k)*(z(i)-x(k))+d(k)
      write(30,*)z(i),s(i)
      end do
      close(30)
      stop
      end program
     

      subroutine calculosplines (n1,xi,xf,l,a,b,c,d,x)
      implicit none
      integer n,i,n1
      real*8 xi,xf,x,a,b,c,d,l,h

      dimension a(0:500),b(0:500),c(0:500),d(0:500),h(0:500)
      dimension x(0:500),l(0:500)
  	  
      n=n1-1
c las h	    
        do i=0,n
          h(i)=abs(xf-xi)/n
        end do
c las d
        do i=0,n
          d(i)=l(i)
        end do 


        call resuelve (n,b,d,h)
clas a        
        do i=1,n
          a(i-1)=(b(i)-b(i-1))/(3.0d0*(x(i)-x(i-1)))
        end do
c las c
        do i=0,n-1
          c(i)=(d(i+1)-d(i))/(x(i+1)-x(i))-(b(i+1)
          
     & +2.0d0*b(i))*(x(i+1)-x(i))/3.0d0
     	end do 
      return 
      end
c en esta subrutina resolveremos el sistema tridiagonal con las ecs de la teoria
      subroutine resuelve (n,b,d,h)
	    implicit none
        integer n,i
        real*8 b,h,d,u,r,ak,bk,ck
      	dimension b(0:500),h(0:500),ak(0:500),d(0:500),bk(0:500)
        dimension ck(0:500),r(0:500),u(0:500)
        
        do i=1,n-1
          r(i)=((3.d0/h(i))*(d(i+1)-d(i)))-((3.d0/h(i-1))*(d(i)-d(i-1)))
	      bk(i)=2.d0*(h(i)+h(i-1))
        end do

        do i=2,n-1
          ak=h(i-1)
        end do
        
        do i=1,n-2
          ck=h(i)
        end do

        do i=2,n-1
          bk(i)=bk(i)-(ck(i-1)*ak(i)/bk(i-1))
          r(i)=r(i)-(r(i-1)*Ak(i)/Bk(i-1))
        end do
	    u(n-1)=r(n-1)/Bk(n-1)
        do i=n-2,1,-1
          u(i)=(r(i)-Ck(i)*u(i+1))/Bk(i)
        end do
        do i=1,n-1
          b(i)=u(i)
        end do
c condiciones naturales
        b(n)=0
        b(0)=0      
      return 
      end
      
        