c Manuel de la Cruz González 70909708H
      program metodoLagrange3
      implicit none
      real*8 x(0:20),f(0:20),k,puntos,width,xf,xi
      integer i
      write(*,*)'Introduzca el punto inicial del intervalo'
      read(*,*)xi
      write(*,*)'Introduzca el punto final del intervalo'
      read(*,*)xf
      write(*,*)'Introduzca el numero de puntos que desea calcular'
      read(*,*)puntos
c Esto pertenece a la parte a) del problema. Debemos sustituir la k antes pedida por 101 k equispaciadas en el intervalo dado.
c      write(*,*)'Introduzca la x de la que quiere calcular su imagen'
c     read(*,*)k

      width=(abs(xf-xi))/(puntos-1)

	  open(21, FILE='prueba.txt')
	  open(27, FILE='puntos.txt')


c Aquí leemos los puntos que nos da el problema.
      open(35,file='datos.txt')
      do i=0,20
        read(35,*)x(i),f(i)
c        print*,x(i),f(i)
      end do
      close(35)
c Aquí comineza el bucle para calcular los 101 puntos de la función interpolada.
      do i=0,puntos-1
        k=xi+width*i
      
        call lagrange (x,f,k)
      end do 
      close(35)
      stop
      end

      subroutine lagrange(r,s,k)
      implicit none
      real*8 r(0:20),s(0:20),a,b,P,k
      integer i,j
      
   	  P=0.0d0
      do i=0,20
c He definido l como a/b porque me daba problemas definiéndolo directamente
        a=1.0d0
        b=1.0d0
       
          do j=0,20
            
            if(i.ne.j)then
              
             a=a*(k-r(j))
             b=b*(r(i)-r(j))
              
            end if
          end do
          
      P=P+((a/b)*s(i))
      end do
c este write lo puese para comprobar que el programa funcionaba bien internamente      
c      write(21,*)a,b,P,s(i) 
c      end do 
c esto es del apartado a).
c      write(*,*)'El valor de P es: ',P
	  write(27,*)k,P
      return
      end