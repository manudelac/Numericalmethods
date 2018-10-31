c Manuel de la Cruz González DNI:70909708H     
      program metodonewton
      implicit none 
      real*8 f,d,beta,x1,x2,dif,maxit
      integer n
      
      write (*,*)'Introduzca el punto inicial'
      read (*,*)x1
      write (*,*) 'Introduzca el valor de la tolerancia'
      read (*,*) beta
      write (*,*)'Introduzca el número máximo de iteraciones'
      read (*,*)maxit
c declaramos la variable n que nos contará el número de iteraciones
      n=0
c asignamos un valor inicial de modo que se cumpla la condición de entrada al bucle      
      dif=abs(x1*3.d0)

      do while (dif.gt.beta)
        n=n+1
c creamos un if por si la derivada en algun punto es 0
        if (d(x1).eq.0.d0) then
          write(*,*)'La derivada en un punto es 0.'
          write(*,*)' No se puede aplicar el metodo de Newton'
        end if
        x2=x1-(f(x1)/d(x1))
        dif=abs(x2-x1)
        x1=x2
c creamos otro if por si el metodo de Newton no converge       
        if (n.gt.maxit) then
          write(*,*)'El metodo no converge'
        end if
      end do

      write(*,*)'La raiz de la funcion es: ',x2
      write(*,*)'El numero de iteraciones es: ',n
      stop
      end program
      
c declaramos la funcion f(x) y la funcion derivada, d(x)


      function f(x)
      implicit none
      real*8 f,x
      f=exp(-x)*sin(1.d0-(x**2))*(1.d0-(x**2))
      return
      end

      function d(x)
      implicit none
      real*8 d,x
      d=-exp(-x)*(sin(1.d0-x**2)*(1-x**2+2*x)+cos(1-x**2)*(2*x-2*x**3))
      return
      end
      