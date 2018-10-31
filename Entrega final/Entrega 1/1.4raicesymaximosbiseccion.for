c Manuel de la Cruz González DNI:70909708H  
      program raicesconmallado
      implicit none 
      real*8 a,b,width,beta,f,c1,c2,d
      integer i,num
c Este programa lo haremos con una subrutina que realice el método de la bisección      
      write(*,*)'introduzca el valor inicial del intervalo'
      read(*,*) a
      write(*,*)' introduzca el valor final del intervalo'
      read(*,*)b
      write(*,*)'introduzca la tolerancia'
      read*,beta
      write(*,*)'introduzca el numero de subintervalos'
      read*,num
      width=(b-a)/num
c el bucle va hasta num-1 porque c2 al final ons produce el intervalo num+1
      do i=0,num-1
        c1=a+width*i
        c2=a+width*(i+1)

        if(f(c1)*f(c2).lt.0) then
          call biseccion(c1,c2,beta)
        end if
        if (d(c1)*d(c2).lt.0) then
          call maximos (c1,c2,beta)
        end if
      end do
      stop
      end program
c creamos la subrutina que realice el metodo de la biseccion y muestre los resultados en pantalla

	  subroutine biseccion (x,y,alfa)
      implicit none
      real*8 x,y,alfa,f,kmenos1,k
      integer j

      j=0
      k=(x+y)/2
      kmenos1=k*3.d0
      do while (abs(k-kmenos1).gt.alfa)
        kmenos1=k
        j=j+1
        if (f(x)*f(k).eq.0.d0)then
          k=k
          else if (f(x)*f(k).lt.0.d0)then
            y=k
            k=(x+k)/2.d0
          else if (f(y)*f(k).lt.0.d0)then
            x=k
            k=(y+k)/2.d0
        end if
      end do
      write(*,*)'el valor de la raiz es: ',k
      write(*,*)'el numero de iteraciones ha sido', j
      return
      end
c creamos una subrutina análoga a la anterior que calcule los máximos y los mínimos
	  subroutine maximos (x,y,alfa)
      implicit none
      real*8 x,y,alfa,f,kmenos1,k
      integer j

      j=0
      k=(x+y)/2
      kmenos1=k*3.d0
      do while (abs(k-kmenos1).gt.alfa)
        kmenos1=k
        j=j+1
        if (f(x)*f(k).eq.0.d0)then
          k=k
          else if (f(x)*f(k).lt.0.d0)then
            y=k
            k=(x+k)/2.d0
          else if (f(y)*f(k).lt.0.d0)then
            x=k
            k=(y+k)/2.d0
        end if
      end do
c este if solo vale para esta función porque podemos ver que los máximos están por encima de 0 y los mínimos por debajo
      if (f(k).gt.0)then
        write(*,*)'el valor del maximo es: ',f(k)
        write(*,*)'el máximo se encuentra en: ',k
        write(*,*)'el numero de iteraciones ha sido', j
      else
        write(*,*)'el valor del mínimo es: ',f(k)
        write(*,*)'el minimo se encuentra en: ',k
        write(*,*)'el numero de iteraciones ha sido', j
      end if
      return
      end
c creamos la funcion para que calcule las imagenes de la ecuacion
     
      function d(x)
      real*8 d,x
      d=-Exp(-x)*((sin(1-x**2))*(1-x**2+2*x)+cos((1-x**2))*(1-x**2)*2*x)
      return 
      end 
c creamos la funcion que calcule las imagenes de la derivada      
   	  function f(x)
      real*8 f,x
      f=exp(-x)*sin(1.d0-(x**2))*(1.d0-(x**2))
      return 
      end 
      