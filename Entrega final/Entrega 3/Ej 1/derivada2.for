c Manuel de la Cruz González 70909708H
      program derivacion2
      implicit none 
      real*8 t,ho,hf,salto,h,ds,da,dr,vr
      real*8 ers,era,err
      integer n,i
c vr es el calculado con mathematica
	  vr=-0.484356249657267
      
c Pedimos lso datos para el problema      
      write(*,*)'Introduzca el valor de t'
      read(*,*)t
      write(*,*)'Introduzca el valor inicial de h'
      read(*,*)ho
      write(*,*)'Introduza el valor final de h'
      read(*,*)hf
      write(*,*)'Introduzca el salto entre valores de h'
      read(*,*)salto

      n=(abs(hf-ho))/salto
c Abrimos un archivo donde escribir los valores y creamos un bucle que rellene el archivo
      open(35,file='derivada.txt')
      
      do i=0,n
        
        h=ho+i*salto

      
c Fórmula simétrica
      call simetrica(t,h,ds,vr,ers)
c Fórmula antisimétrica      
      call antisimetrica(t,h,da,vr,era)
c Extrapolación de Richardson      
      call richardson(t,h,dr,vr,err)
c Escribimos en un fichero las derivadas      
      write(35,99)h,ds,da,dr,ers,era,err
c Mostramos en pantalla los valores      
      write(*,99)h,ds,da,dr,ers,era,err
      
      end do
c Formato del fichero que hemos abierto      
99 	  format(f10.3,1x,3(f20.15,3x),3(e10.4,3x))   
      close(35)
      stop
      end

c Subrutina simétrica     
      subroutine simetrica(t,h,ds,vr,ers)
      implicit none
      real*8 t,h,v,ds,ers,vr
      

      ds=(v(t+h)-v(t-h))/(2*h)
c      write(*,*)'el valor de la derivada por el met sim en t es: ',ds
c Introducimos en todas las subrutinas su error relativo correspondiente (los er)
	  ers=((ds-vr)/vr)*100.d0
      return
      end
    
c Subrutina antisimétrica
      subroutine antisimetrica(t,h,da,vr,era)
      implicit none
      real*8 t,h,v,da,era,vr
      
	  da=(v(t+h)-v(t))/h
	  era=((da-vr)/vr)*100.d0
c      write(*,*)'el valor de la derivada por met anti en t es: ',danti

      return 
      end
c c Subrutina extrapolación de Richardson
      subroutine richardson(t,h,dr,vr,err)
      implicit none
      real*8 t,h,v,dr,err,vr

      dr=(-1.0d0/(6.0d0*h))*(v(t+h)-8.0d0*v(t+(h/2.0d0))-v(t-h)+
     & 8.0d0*v(t-(h/2)))
      err=((dr-vr)/vr)*100.d0
     
c      write(*,*)'el valor de la derivada por interp rich en t es: ',rich
     
	  return
      end
c FUNCIÓN PROBLEMA
      function v(t)
      implicit none
      real*8 t,v

      v=(dexp(-t))*dcosh((1.0d0-t)**2)

      return 
      end