program calorCN
    implicit none

    ! parametres fisics
    real(8), parameter :: L = 0.02d0
    real(8), parameter :: k = 0.56d0
    real(8), parameter :: rho = 1081.d0
    real(8), parameter :: cv = 3686.d0
    real(8), parameter :: Pext = 9.44d5
    real(8), parameter :: Tref = 309.65d0
    real(8), parameter :: deltaT = 43.5d0
    real(8), parameter :: Ttarget = 353.15d0

    ! escalar adimensional
    real(8), parameter :: t0 = (L*L*rho*cv)/k
    real(8), parameter :: C = Pext*L*L/(k*deltaT)

    ! parametres numerics
    integer, parameter :: Nx = 200
    real(8), parameter :: dx = 1.d0/(Nx-1)
    real(8), parameter :: dt = 1.d-5
    integer, parameter :: max_steps = 2000000

    real(8) :: r, Tmax, t_dim
    integer :: i, step

    ! arrays del metode tridiagonal
    real(8), dimension(Nx) :: theta, theta_new
    real(8), dimension(Nx) :: aa, bb, cc, dd

    ! inicialitzacio
    r = dt/(2.d0*dx*dx)

    ! condicio inicial: tot zero
    theta = 0.d0

    ! inicialitzar matriu tridiagonal amb identitat
    do i = 1, Nx
        aa(i) = 0.d0
        bb(i) = 1.d0
        cc(i) = 0.d0
    end do

    ! assignar coeficients Crank-Nicolson als punts interiors
    do i = 2, Nx-1
        aa(i) = -r
        bb(i) = 1.d0 + 2.d0*r
        cc(i) = -r
    end do

    ! bucle temporal
    do step = 1, max_steps

        ! condicions de frontera al vector costat dret
        dd(1) = 0.d0
        dd(Nx) = 0.d0

        ! construir vector costat dret segons Crank-Nicolson
        do i = 2, Nx-1
            dd(i) = r*theta(i-1) + (1.d0 - 2.d0*r)*theta(i) + r*theta(i+1) + C*dt
        end do

        ! resoldre sistema tridiagonal amb Thomas
        call thomas(aa, bb, cc, dd, theta_new, Nx)
        theta = theta_new

        ! calcular temperatura maxima real
        Tmax = maxval(theta)*deltaT + Tref

        ! sortir si s'ha arribat a la temperatura objectiu
        if (Tmax >= Ttarget) exit
    end do

    ! calcular temps dimensional
    t_dim = step * dt * t0

    ! impressio dels resultats
    print *, "resultats"
    print *, "temps per arribar a 80 graus C:", t_dim, " s"
    print *, "temperatura maxima final:", Tmax, " K"
    print *, "posicio (m)   temperatura (K)"

    do i = 1, Nx
        print *, (i-1)*L/(Nx-1), theta(i)*deltaT + Tref
    end do

end program calorCN

! rutina thomas per resoldre sistemes tridiagonals
subroutine thomas(a, b, c, d, x, n)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(in) :: a, b, c
    real(8), dimension(n), intent(in) :: d
    real(8), dimension(n), intent(out) :: x

    real(8), dimension(n) :: cp, dp
    real(8) :: m
    integer :: i

    ! convertim el sistema en un triangular superior
    ! amb aquestes transformacions

    cp(1) = c(1) / b(1)
    dp(1) = d(1) / b(1)

    do i = 2, n
        m = b(i) - a(i)*cp(i-1)
        cp(i) = c(i) / m
        dp(i) = (d(i) - a(i)*dp(i-1)) / m
    end do

    ! solucionem a partir de x_n

    x(n) = dp(n)
    do i = n-1, 1, -1
        x(i) = dp(i) - cp(i)*x(i+1)
    end do

end subroutine thomas