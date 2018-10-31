c Manuel de la Cruz González 70909708H
	  program sistemamatrices
      implicit none
      real*8 a1,b1,x2(1:100),a2,b2,x1 (1:100),det1,det2
c   a1,b1,a2 y b2 son las matrices. x2 y x1 los vectores de incognitas.   
      integer i,j,c,n
c i y j son indices de bucles
c c es el numero que selecciona el case
c n es la dimension de la matriz      
      dimension a1(100,100),b1(1:100),a2(100,100),b2(1:100)
      
      write(*,*)'Introduzca la dimension de la matriz'
      read*, n   	
      n=4
      open(35,file='a1.txt')
      open(36 ,file='b1.txt')
c Aquí leemos la matriz      
      do i=1,n
        read(35,*)(a1(i,j),j=1,n)
        read(36,*)b1(i)
      end do
      close (35)
      close (36)

     
      open(37,file='a2.txt')
      open(38 ,file='b2.txt')
      do i=1,n
        read(37,*)(a2(i,j),j=1,n)
        read(38,*)b2(i)
      end do
      close (37)
      close (38)
      
c Inicializamos a 0s los vectores
      do i=1,n
        x2(i)=0.d0
        x1(i)=0.d0
      end do        
c Para comprobar que ha leido bieno la matriz. Poner despues de a el número que corresponde.
c Hay que añadir el numero a a y a b
c      do i=1,n
c       	write(*,'(10f8.3)') (a(i,j), j=1,n)
c      end do
c      print*
c      write(*,'(10f8.3)')(b(i),i=1,n)

c Abrimos el case para que seleccione si es triangular superior o inferior.
	  write(*,*)'Si la matriz es triangular superior introduzca 0'  
      write(*,*)'Si la matriz es triangular inferior introduzca 1'  
	  read*,c
      select case(c)
        case(0)
              call sdetras (a2,b2,n,x2)
          do i=1,n
            write(*,*) x2(i)
          end do 
            call determ (a2,n,det1)
            write(*,*)'El determinante de la matriz triangular sup es '
            write(*,*) det1 
        case(1)
            call sdelante(a1,b1,n,x1)
            do i=1,n
              write(*,*) x1(i)
            end do 
            call determ (a1,n,det2)
            write(*,*)'El determinante de la matriz triangular inf es '
            write(*,*) det2 
         
        case default
          print*,'Error del sistema'
      end select
	
 

      stop
      end
c sustitución hacia atrás se hace en matrices triangulares superiores
      subroutine sdetras (a,b,k,x2)
	  implicit none
      real*8 sum,a(100,100),b(1:100),x2(1:100)
      integer k,i,j
      sum=0
      x2(k)=b(k)/a(k,k)
      x2(k-1)=(b(k-1)-a(k-1,k)*x2(k))/(a(k-1,k-1))
      
      do i=k-2,1,-1
        
        do j=i+1,k
	    	sum=sum+a(i,j)*x2(j)            
        end do    	
      x2(i)=(b(i)-sum)/a(i,i)  
      end do
      
      return 
      end

c sustitución hacia delante se hace en matrices triangulares inferiores      
      subroutine sdelante (a,b,k,x1)
	  implicit none
      real*8 sum,a(100,100),b(1:100),x1(1:100)
      integer k,i,j
      sum=0
      x1(1)=(b(1)/a(1,1))
      x1(2)=(b(2)-a(2,1)*x1(1))/(a(2,2))
      do i=3,k
        do j=1,i-1
          sum=sum+a(i,j)*x1(j)
        end do  
      x1(i)=(b(i)-sum)/a(i,i)   	
        
      end do
      
      return 
      end

c Subrutina para calcular los determinantes  

      subroutine determ (a,k,det)
	  implicit none
      real*8 a(100,100),det
      integer k,i,j
      det=1.d0
      do i=1,k
        do j=1,k
          if (i.eq.j) then
            det=det*a(i,j)
          end if
        end do
      end do
      return
      end      
         