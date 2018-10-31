c Manuel de la Cruz González DNI:70909708H  
      PROGRAM biseccion
      implicit none
      real*8 a,b,c,beta,cmenos1,f
      integer n

      write (*,*)'Introduzca el intervalo separado por comas'
      read (*,*)a,b
      write (*,*) 'Introduzca el valor de la tolerancia'
      read (*,*) beta
c declaramos una variable n que nos contará las iteracciones
      n=0
      c=(a+b)/2.d0
      cmenos1=c+1
c la variable cmenos1 tomará el valor de c una vuelta antes en el bucle
      do while (abs(c-cmenos1).gt.beta)
        cmenos1=c
        n=n+1
c escribimos un if para que en el caso de que a o b sean raíces salga del bucle
        if (f(a)*f(b).eq.0.d0)then
          write(*,*)'Introduzca un intervalo nuevo'
          else if (f(a)*f(c).lt.0.d0)then
            b=c
            c=(a+c)/2.d0
          else if (f(b)*f(c).lt.0.d0)then
            a=c
            c=(b+c)/2.d0
        end if
      end do 
      write(*,*) 'El valor de la raíz para su tolencia es c',c
      write (*,*)'El numero de iteracciones ha sido: ',n
      stop
      end program
c creamos una función para calcular f(_)
	  function f(x)
      real*8 f,x
      f=exp(-x)*sin(1.d0-(x**2))*(1.d0-(x**2))
      return 
      end
            
          














     

      