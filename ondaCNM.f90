program onda
! Usando el metodo de Crank-Nocolson se va a aproximar la ecuacion de onda: ∂^2u/∂t^2=∂^2u/∂x^2
! Con las condiciones de frontera:  0<x<1 ; h=0.1 ; k=0.1
! u(0,t)=u(1,t)=0 ; u(x,0)=sen(2PiX) ; ∂u/∂t(x,0)=2Pisen(2Pix)
    real h,k
    real xl,xu,tl,tu
    real x,t
    integer contadorx,contadort
    real, dimension (:,:),ALLOCATABLE :: m !REAL, DIMENSION(:), ALLOCATABLE

    !paso de integracion  h:espacio k:tiempo
    h=0.1
    k=0.1
    !rango de integracion
    xl=0.0
    xu=1.0
    tl=0.0
    tu=10.0
    !inicializacion de las coordenadas
    x=xl
    t=tl
    !contadores
    contadorx=0
    contadort=0
    !Creer un archivo para guardar los datos
    open(22,file="archivoOnda.txt",status='new')

    do !contador en x
        x=x+h
        if(x>xu) then
            x=xl
            exit
        end if
        contadorx=contadorx+1
    end do

    do !contador en t
        t=t+k
        if(t>tu) then
            t=tl
            exit
        end if
        contadort=contadort+1
    end do
    !Los contadores se utilizan para indicar las dimenciones para la matriz
    allocate( m(0:contadorx,0:contadort))


    m(0,0)=u(0.0)
    m(contadorx,0)=u(xu)

    do i=0 , contadort
        m(0,i)=0.0
        m(1,i)=0.0
    end do

    do i=1, contadorx-1
        m(i,1)=ui1(i*h,k)
    end do

    do j=1 , contadort-1
        do i=1 , contadorx-1
            m(i,j+1)=(1.0/3.0)*(3.0*m(i+1,j)-m(i,j-1)+m(i-1,j)-m(i+1,j-1)+m(i-1,j+1))
        end do
    end do

    do j=0 , contadort
        do i=0 , contadorx
            x=i*h
            t=j*k
            print*, t,x,m(i,j)
            write(22,'(F12.7,F12.7,F15.11)') t,x,m(i,j) !guardar datos
        end do
    end do
end program


    real function u (x) result(resultado)!funcion que retorna alguna de las condiciones u(0,t)=u(1,t)=0 ; u(x,0)=sen(2PiX)
        real x
        real Pi
        !constante π
        Pi=3.1415927
        resultado=sin(2.0*Pi*x)
    end function

    real function ui1(x,k) result(resultado1) !funcion con la forma u(xi,1)=kg(x)+u(xi,0) se utiliza la condicion de frontera ¶u/¶t(x,0)=2Pisen(2Pix)
        real g
        real Pi
        real k
        !constante π
        Pi=3.1415927
        g=2.0*Pi*sin(2.0*Pi*x)
        resultado1=k*g+u(x) ! se puede usar una aproimacion mas apropiada
    return
    end function
