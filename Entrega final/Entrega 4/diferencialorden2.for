c Manuel de la Cruz González      70909708H
	  program diferencialorden2
      implicit none
      real*8 f2,f1,ti,tf,h,q0,i0,t,i,q,k11,k12,k13,k14,k21,k22,k23,k24
      real*8 r,l,c
      integer n,j
      dimension i(0:200),q(0:200),t(0:200)
c Pedimos los valores para el problema
      write(*,*)'Valor de R'
      read(*,*)r
      write(*,*)'Valor de L'
      read(*,*)l
      write(*,*)'Valor de C'
      read(*,*)c
	
      write(*,*)'Inicio del intervalo de tiempo'
      read(*,*)ti
      write(*,*)'Final del intervalo de tiempo'
      read(*,*)tf
	  
      write(*,*)'Numero de puntos'
      read(*,*)n
      
      write(*,*)'Condiciones iniciales de i (t=0)'
      read(*,*)i0
	  
      write(*,*)'Condiciones iniciales de q (t=0)'
      read(*,*)q0
      
      h=(tf-ti)/(n-1)
      open(22,file='Q(t).txt')
      open(23,file='I(t).txt')
      
	  i(0)=i0
      q(0)=q0
      t(0)=ti
c Bucle que calcula los R4K en diferentes valores de t. Calculamos un runge kutta para cada funcion.
c F1 será la derivada segunda de q respecto de t y f2 dq/dt, la intensidad

      do j=1,n-1
        t(j)=ti+j*h
        k11=f1(i(j-1),q(j-1),r,l,c)
        k21=f2(i(j-1))
        k12=f1(i(j-1)+(h/2.d0)*k11,q(j-1)+(h*k21/2.d0),r,l,c)
        k22=f2(i(j-1)+(h/2.d0)*k11)
        k13=f1(i(j-1)+(h/2.d0)*k12,q(j-1)+(h*k22/2.d0),l,c,r)
        k23=f2(i(j-1)+(h/2.d0)*k12)
        k14=f1(i(j-1)+h*k13,q(j-1)+h*k23,r,l,c)
        k24=f2(i(j-1)+h*k13)
        
        i(j)=i(j-1)+(h/6.d0)*(k11+2.d0*k12+2.d0*k13+k14)
        q(j)=q(j-1)+(h/6.d0)*(k21+2.d0*k22+2.d0*k23+k24)
      end do
c Aquí escribimos los datos obtenidos en cada t.
      do j=0,n-1
        write(22,*) t(j),q(j)
      end do
      
      do j=0,n-1
        write(23,*) t(j),i(j)
      end do
      
      close(22)
      close(23)
      
      stop
      end

      
c      FUNCIONES

      function f1(i,q,r,l,c)
      implicit none
      real*8 f1,r,q,l,c,i


      f1=-(r*c*i+q)/(l*c)

      return
      end


      function f2(i)
      implicit none
      real*8 i,f2

      f2=i

      return
      end

      