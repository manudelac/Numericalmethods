	  program gaussseidelsinpivot
      implicit none
      real*8 A,b,x,e
      integer i,n,j
      dimension  A(100,100), b(100), x(100)

	  write(*,*)'Orden de la matriz: '
      read*,n	
      print*
      write(*,*)'Introduzca la precisión'
      read*,e

c Abrimos los ficheros con A y b     
      open(30,file='A.txt')
      open(31,file='b.txt')
      open(32,file='x.txt')
      
c Leemos en Fortran A y b     
      do i=1,n
         read(30,*)(A(i,j),j=1,n)
      end do
      do i=1,n
         read(31,*) b(i)
      end do
      close(30)
      close(31)
c Escribimos las condiciones iniciales
      x(1)=0.01d0
      x(2)=0.01d0
      x(3)=0.01d0
      x(4)=0.01d0
      x(5)=0.01d0
c LLamada a la subrutina   
      call gaussseidel(a,b,x,n,e)
c Escribimos x en un archivo externo     
      do i=1,n
        write(32,*) x(i)
      end do

      stop
      end
      
c SUBRUTINA

	  subroutine gaussseidel(a,b,x,n,e)
      implicit none
      real*8 A,b,x,e,converg,sum1,sum2,xantes,paso
      integer i,j,n,k
      dimension A(100,100), b(100), x(100), xantes(100)
c Tenemos que crear una variable que cuente las iteraciones , será k.
	  k=0      
c converg tiene que ser inicializada un numero grande para que entre en el bucle.
	  converg=e+1.d0
c como este método tiene como criterio de stop uno de convergencia, usamos un bucle do while con la cond de convergencia
      do while(converg.ge.e)
        converg=0.d0
c sumamos una vuelta por iteracion        
        k=k+1
c aqui empiezan los do
        do i=1,n
          do j=1,n
          xantes(j)=x(j)
        end do
        sum1=0.d0  
        do j=1,i-1
          sum1=sum1+A(i,j)*x(j)
        end do
        sum2=0.d0
        do j=i+1,n
          sum2=sum2+A(i,j)*xantes(j)
        end do
        x(i)=1.d0/A(i,i)*(b(i)-sum1-sum2)    
        paso=abs(x(i)-xantes(i))
c escribimos una variable que nos guarde la dif en valor absoluto en xi de k y xi de k+1        
          paso= abs(x(i)-xantes(i))

	      if (paso.gt.converg) then
      	  converg=paso
          end if
        end do
      end do
      write(*,*)'El numero de iteraciones ha sido',k
      return
      end        



    