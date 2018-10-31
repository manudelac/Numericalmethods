c Manuel de la Cruz González 70909708H
	  program ludoolittle
      implicit none
      real*8 A,L,U,det
      integer i,j,n
      dimension A(100,100),L(100,100),U(100,100)


	  write(*,*)'Introduzca el orden de la matriz'
      read*,n
c creamos archivos en los que obtener las matrices que obtenemos

      open(33,file='A.txt')
      open(34,file='matrizL.txt')
      open(35,file='matrizU.txt')      
c Leemos la matriz desde el fichero de texto
  	  do i=1,n
         read(33,*)(A(i,j),j=1,n)
      end do    
c Inicializamos a cero las matrices L y U
	  do i=1, n
        do j=1, n
          L(i,j)=0.d0
          U(i,j)=0.d0
        end do
      end do    
c Llamamos a la subrutina LU para que opere con A.
	  call LU(A,L,U,n)
      
c Escribimos en los archivos las matrices que obtenemos 
	  do i=1,n
         write(34,999)(L(i,j),j=1,n)
         write(35,999)(U(i,j),j=1,n)
999      format(5 d18.8)
      end do      
     
   
c det A es el producto de detU y det L     
      det=1.d0
      do i=1,n
         det=det*L(i,i)*U(i,i)
      end do  
      write(*,*)'El determinante de A es: ',det
              
      stop
      end
c Subrutina LU

	  subroutine LU(A,L,U,n)
      implicit none
      real*8 A,L,U,sum1,sum2
      integer n,i,j,k,t
      dimension A(100,100),L(100,100),U(100,100)

       
      do i=1,n
        L(i,i)=1
        U(1,i)=A(1,i)
      end do  

	  do k=2,n
      	do i=1,k-1
        	L(k,i)=0.d0
        end do
      end do   

      do k=2,n
        do i=1,k-1
          sum1=0.d0
          do j=1,i-1
            sum1=sum1+L(k,j)*U(j,i)
          end do
        L(k,i)=(A(k,i)-sum1)/U(i,i)
        end do   

         do t=k,n 
           sum2=0.d0         
            do j=1,t-1
              sum2=sum2+L(k,j)*U(j,t)
            end do  
            U(k,t)=A(k,t)-sum2
         end do  
       end do

      return
      end



     