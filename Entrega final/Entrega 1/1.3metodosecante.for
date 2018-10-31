c Manuel de la Cruz González DNI:70909708H  
      program metodosecante
      implicit none
      real*8 x1k,xk,xk1,f,beta,dif,maxit
      integer n

      write(*,*)'introduzca xk-1'
      read(*,*)x1k
      write(*,*)'introduzca xk'
      read(*,*)xk
      write (*,*) 'Introduzca el valor de la tolerancia'
      read (*,*) beta
      write (*,*) 'Numero maximo de iteraciones'
      read (*,*) maxit
      
      
c declaramos una variable que cuente el numero de iteraciones
	  n=0
c asignamos un valor inicial de modo que se cumpla la condición de entrada al bucle     
      dif=abs(beta*3.d0)
c el bucle calculará xk1(x_k+1)y      

      do while (dif.gt.beta)
 
        n=n+1
        xk1=xk-f(xk)*((xk-x1k)/(f(xk)-f(x1k)))
        dif=abs(xk1-xk)
        x1k=xk
        xk=xk1
c introducimos un bucle if por si el metodo no converge       
        if (n.gt.maxit) then
            write(*,*)'el método de la secante no converge'
        end if
      end do

      write(*,*)'La raiz encontrada de la funcion es: ',xk1
      write(*,*)'el numero de iteracciones ha sido:  ',n

      stop
      end program
c creamos una funcion que calcule las imagenes
	  function f(x)
      implicit none
      real*8 f,x
      f=exp(-x)*sin(1.d0-(x**2))*(1.d0-(x**2))
      return
      end

      
        
      

      