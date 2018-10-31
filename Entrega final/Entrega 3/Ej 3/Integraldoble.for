c Manuel de la Cruz González 70909708H
      program integraldoble
      implicit none
      real*8 ti,tf,a(0:20),b(0:20),k,intfinal,erint,vr
      integer i
c Valor calculado en mathematica   
      vr=2.63502222757
c Archivo con los pesos y ceros de legendre para 10 puntos
      open(22,file='legendre.txt')
c Leemos el archivo con los pesos y los ceros y los metemos en vectores        
      do i=0,9   
      read(22,*) a(i),b(i)
      end do
      write(*,*)'Extremo inicial del intervalo de T'
      read*, ti
      write(*,*)'Extremo inicial del intervalo de T'
      read*, tf
c Llamamos a una subrutina, que nos dara el valor de la integral en los intervalos dados 
	  call integraltemp(a,b,ti,tf,k)
      intfinal=0.5*(tf-ti)*k
      print*, 'el valor de la integral es', intfinal
      erint=Abs(intfinal-vr)/vr
      print*, 'el error relativo es:', erint
      stop
      end
c Aquí calculamos la t como una constante y se la pasamos a otra subrutina para que con esa información haga la integral doble. 

      subroutine integraltemp(a,b,ti,tf,k)
      implicit none
      integer i
      real*8 ti,tf,a(0:20),b(0:20),p,k,tc,dif
      dif=tf-ti
      k=0.d0
      do i=0,9
        tc=0.5*((tf+ti)+(dif*b(i)))
        call integralp(b,a,tc,p)
        k=k+a(i)*p
      end do
      return
      end
     
c Aquí en cada ciclo del bucle se calculará p mediante legendre de 10 puntos y volverá a la subrutina anterior

      subroutine integralp(b,a,tc,p)
      implicit none
      real*8 ti,tf,a(0:20),b(0:20),h,tc,p,f1,f2,p1,k
      integer i
      ti=f1(tc)
      tf=f2(tc)
      p=0.d0
      do i=0,9
        P1=0.5*((tf+ti)+((tf-ti)*b(i)))
        p=p+a(i)*h(tc,p1)
      end do
      p=p*0.5*(tf-ti)
      return
      end
      
c FUNCIONES QUE RELACIONAN T Y P     
      function f1(T)
      real*8 f1,t
      f1=0.d0
      return
      end

      function f2(T)
      implicit none
      real*8 f2,T
      f2=1.d0+T
      return
      end
      
c FUNCION
      function h(T,P)
      real*8 h,T,P
      h=((T**2)+P**2)*exp(-T*P)
      return
      end