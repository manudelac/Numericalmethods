c Manuel de la Cruz González   70909708H
	  program ecuaciondiferencial
      implicit none
      real*8 xi,xf,x,y,width,f,K1,K2,K3,K4
      integer np,k,i

      write(*,*)'Introduzca el extremo inicial del intervalo'
      read(*,*)xi
	  
      write(*,*)'Introduzca el extremo final del intervalo'
      read(*,*)xf
	  
      write(*,*)'Introduzca el numero de pares'
      read(*,*)np
	  
	  k=np-1
      width=(abs(xf-xi))/k

c Abrimos un archivo para cada método donde escribiremos los datos obtenidos

	  open(12,file='Taylorprimer.txt')
      open(13,file='Heun.txt')
      open(14,file='Eulermodificado.txt')
      open(15,file='Rugekutta4.txt')
c Método de Taylor de primer orden 
	  y=0.d0
      x=0.d0
      write(12,*)x,y
      do i=1,k
        y=y+width*f(x,y)
        x=i*width
      write(12,*)x,y
      end do
      close(12)
      
c Método de Heun
      y=0.d0
      x=0.d0
      write(13,*)x,y
      do i=1,k
        y=y+(width/2.d0)*(f(x,y)+f(x+width,y+width*f(x,y)))
        x=i*width
      write(13,*)x,y
      end do
      close(13)
      

c Método de euler modificado
      y=0.d0
      x=0.d0
      write(14,*)x,y
      do i=1,k
        y=y+width*f(x+(width/2.d0),y*(width/2.d0)*f(x,y))
        x=i*width
      write(14,*)x,y
      end do
      close(14)

c Método Runge-Kutta de orden 4
	  y=0.d0
      x=0.d0
      write(15,*)x,y 
      do i=1,k
        k1=f(x,y)
        k2=f(x+width/2.d0,y+(width/2.d0)*k1)
        k3=f(x+width/2.d0,y+(width/2.d0)*k2)
        k4=f(x+width,y+width*k3)
        y=y+(width/6.d0)*(k1+2.d0*k2+2.d0*k3+k4)
        x=i*width
        write(15,*)x,y
      end do       
      stop
      end      

c FUNCION
	  function f(x,y)
      real*8 y,x,f
      f=y-exp(x)*cos(x)*sin(x) 
      return
      end