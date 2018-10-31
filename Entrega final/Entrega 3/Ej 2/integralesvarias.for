c Manuel de la Cruz González 70909708H     
      program integrales
      implicit none i
      real*8 xi,xf,pm,trap,simp,new,gl2,gl10,vr,erpm,ertr,ersp,ernw
      real*8 ergl2,ergl10,trapcomp,ertc,simpcomp,ersc,n,n2

      write(*,*)'Introduzca el intervalo separados por comas :'
      read (*,*)xi,xf
c Aquí pedimos el número de puntos para realizar las reglas compuestas.
      write(*,*)'Numero de puntos para el trapecio compuesto :'
	  read(*,*) n
      write(*,*)'Numero de puntos para simpson compuesto :'
	  read(*,*) n2

      xi=0.0d0
      xf=1.8d0
c vr es el valor de la integral calculado en mathematica, lo utilizaremos para medir el error relativo     
      vr=0.930797305275
C Llamamos a las subrutinas que calculan por cada método
C Punto medio
	  call puntomedio(xi,xf,pm,erpm,vr)
      print*
      write(*,*)'Punto medio'
      write(*,*)'Valor',pm
      write(*,*)'Error relativo:',erpm
c Trapecio simple      
      call trapecio(xi,xf,trap,ertr,vr)
      print*
      write(*,*)'Trapecio simple'
      write(*,*)'Valor',trap
      write(*,*)'Error relativo:',ertr
c Regla de Simpson     
      call simpson (xi,xf,simp,ersp,vr)
      print*
      write(*,*)'Regla de Simpson'
      write(*,*)'Valor',simp
      write(*,*)'Error relativo:',ersp
c Newton-Cotes      
      call newtoncotes(xi,xf,new,ernw,vr)
      print*
      write(*,*)'Newton-Cotes'
      write(*,*)'Valor',new
      write(*,*)'Error relativo:',ernw
c Gauss-legendre dos puntos
      call gausslegendre2(xi,xf,gl2,ergl2,vr)
      print*
      write(*,*)'Gauss-Legendre 2 puntos'
      write(*,*)'Valor',gl2
      write(*,*)'Error relativo:',ergl2
c Gauss-legendre 10 puntos
      call gausslegendre10(xi,xf,gl10,ergl10,vr)
      print*
      write(*,*)'Gauss-Legendre 10 puntos'
      write(*,*)'Valor',gl10
      write(*,*)'Error relativo:',ergl10
c Trapecio compuesto
  	  call trapeciocompuesto(xi,xf,trapcomp,n,ertc,vr)
      print*
      write(*,*)'Trapecio compuesto'
      write(*,*)'Valor',trapcomp
      write(*,*)'Error relativo:',ertc
c Simpson compuesta
	  call simpsoncomp(xi,xf,simpcomp,n2,ersc,vr)
      print*
      write(*,*)'Simpson compuesta'
      write(*,*)'Valor',simpcomp
      write(*,*)'Error relativo:',ersc
	  stop
      end
      
c SUBRUTINAS
      subroutine puntomedio(xi,xf,pm,erpm,vr)
c n=0      
      implicit none
      real*8 xf,xi,pm,v,vr,erpm
      pm=(xf-xi)*v((xi+xf)/2.0d0)
      erpm=abs((vr-pm)/vr)
      return
      end

	  subroutine trapecio(xi,xf,trap,ertr,vr)
c n=1      
      implicit none
      real*8 xf,xi,trap,h,v,ertr,vr
      h=(xf-xi)
      trap=(h/2.0d0)*(v(xi)+v(xf))
      ertr=abs((vr-trap)/vr)
      return
      end

      subroutine simpson(xi,xf,simp,ersp,vr)
c n=2      
      implicit none
      real*8 xf,xi,v,simp,vr,ersp
      
      simp=((xf-xi)/6.0d0)*(v(xi)+(4.0d0*v((xi+xf)/2.0d0)+v(xf)))
      ersp=abs((vr-simp)/vr)
      return
      end

      subroutine newtoncotes(xi,xf,new,ernw,vr)
      implicit none
      real*8 xf,xi,v,new,h,vr,ernw
      h=(xf-xi)/4.d0
      new=0.14*v(xi)+0.64*v(xi+h)+0.24*v(xi+2*h)+0.64*v(xi+3*h)+
     & 0.14*v(xf)
      
      return
      end
c Para la integral mediante gauss-legendre con 2 puntos lo escribimos directamente de la fórmula
      subroutine gausslegendre2(xi,xf,gl2,ergl2,vr)
      implicit none
      real*8 h,xi,xf,x0,x1,gl2,alfa0,alfa1,v,vr,ergl2
      h=xf-xi
      alfa0=1.0d0
      alfa1=alfa0
      x0=0.577350269189626d0
      x1=-0.577350269189626d0
      gl2=(h/2.0d0)*(alfa0*v((xi+xf+(h*x0))/2.d0)+alfa1*v((xi+xf+(h*x1))
     &/2.d0)) 
      ergl2=abs((vr-gl2)/vr)
      return
      end
c Para la integral mediante gauss legendre con 10 puntos creamos un archivo con los puntos y le pedimos que lo lea

      subroutine gausslegendre10(xi,xf,gl10,ergl10,vr)
      integer i
      real*8 alfa(0:10),equis(0:10),sum,v,gl10,xi,xf,vr,ergl10
      sum=0
      do i=0,9
       open(38,file='legendre.txt')
       read(38,*) alfa(i),equis(i)
      sum=sum+alfa(i)*v((xf+xi+(xf-xi)*equis(i))/2)
      end do
      gl10=(1.0d0/2.0d0)*(xf-xi)*sum
      ergl10=abs((vr-gl10)/vr)
      return
      end
c REGLAS COMPUESTAS
c Regla del trapecio compuesto.  
	  subroutine trapeciocompuesto(xi,xf,trapcomp,n,ertc,vr)
      implicit none
      integer i
      real*8 xi,xf,trapcomp,n,ertc,vr,h,k,sum,wi,v
      h=xf-xi
      wi=xi
      k=(xf-xi)/(n-1)
      sum=0
      do i=1,n-1
      sum=sum+v(wi)
      wi=xi+k*i
      end do
      trapcomp=(h/(n-1))*(((v(xi)+v(xf))/2.d0)+sum)
      ertc= abs((vr-trapcomp)/vr)
      return 
      end
c Regla de Simpson compuesta     
      subroutine simpsoncomp(xi,xf,simpcomp,n2,ersc,vr)
      implicit none
      integer i
      real*8 xi,xf,h,sum1,sum2,k,n2,wi,wi2,simpcomp,v,ersc,vr
      h=xf-xi
      sum1=0
      sum2=0
      k=(xf-xi)/(n2-1)
      do i=1,(n2-1)/2
        wi=xi+(2*i-1)*k
      sum1=sum1+v(wi)
      end do
      do i=1,((n2-1)/2)-1
        wi2=xi+2*i*k
      sum2=sum2+v(wi2)
      end do
      simpcomp=(h/(3*n2))*(v(xi)+v(xf)+4*sum1+2*sum2)
      ersc=abs((vr-simpcomp)/vr)
      return
      end
      
      
      



c FUNCIÓN
      function v(t)
      implicit none
      real*8 t,v

      v=(dexp(-t))*dcosh((1.0d0-t)**2)

      return 
      end