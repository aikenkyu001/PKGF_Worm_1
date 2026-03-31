module pkgf_math_utils
    use pkgf_constants
    implicit none
contains
    subroutine invert_matrix(A, Ainv)
        real(DP), intent(in) :: A(DIM, DIM)
        real(DP), intent(out) :: Ainv(DIM, DIM)
        real(DP) :: tA(DIM, DIM), f, pivot
        integer :: i, j, k, p
        tA = A
        Ainv = 0.0d0
        do i = 1, DIM
            Ainv(i, i) = 1.0d0
        end do
        do i = 1, DIM
            p = i
            pivot = abs(tA(i, i))
            do j = i + 1, DIM
                if (abs(tA(j, i)) > pivot) then
                    p = j
                    pivot = abs(tA(j, i))
                end if
            end do
            if (p /= i) then
                do k = 1, DIM
                    f = tA(i, k); tA(i, k) = tA(p, k); tA(p, k) = f
                    f = Ainv(i, k); Ainv(i, k) = Ainv(p, k); Ainv(p, k) = f
                end do
            end if
            f = tA(i, i)
            ! 【浄化】安定化のための微小加算を全廃。特異点は特異点として扱う。
            tA(i, :) = tA(i, :) / f
            Ainv(i, :) = Ainv(i, :) / f
            do j = 1, DIM
                if (i /= j) then
                    f = tA(j, i)
                    tA(j, :) = tA(j, :) - f * tA(i, :)
                    Ainv(j, :) = Ainv(j, :) - f * Ainv(i, :)
                end if
            end do
        end do
    end subroutine invert_matrix

    function determinant(A) result(det)
        real(DP), intent(in) :: A(DIM, DIM)
        real(DP) :: det, tA(DIM, DIM), f
        integer :: i, j
        tA = A
        det = 1.0d0
        do i = 1, DIM
            if (abs(tA(i, i)) < 1.0d-30) then
                det = 0.0d0
                return
            end if
            det = det * tA(i, i)
            do j = i + 1, DIM
                f = tA(j, i) / tA(i, i)
                tA(j, i:DIM) = tA(j, i:DIM) - f * tA(i, i:DIM)
            end do
        end do
    end function determinant

    subroutine matrix_exp(A, expA)
        real(DP), intent(in) :: A(DIM, DIM)
        real(DP), intent(out) :: expA(DIM, DIM)
        real(DP) :: A2(DIM, DIM), A4(DIM, DIM), A6(DIM, DIM)
        real(DP) :: U(DIM, DIM), V(DIM, DIM), num(DIM, DIM), den(DIM, DIM), id(DIM, DIM), den_inv(DIM, DIM)
        integer :: i
        id = 0.0d0
        do i = 1, DIM; id(i, i) = 1.0d0; end do
        A2 = matmul(A, A); A4 = matmul(A2, A2); A6 = matmul(A4, A2)
        U = A * matmul(A6 + 15120.0d0 * A4 + 604800.0d0 * A2 + 17643225600.0d0 * id, id) * (1.0d0/30240.0d0)
        V = 30240.0d0 * A6 + 3326400.0d0 * A4 + 86486400.0d0 * A2 + 17643225600.0d0 * id
        num = V + U; den = V - U
        call invert_matrix(den, den_inv)
        expA = matmul(den_inv, num)
    end subroutine matrix_exp

    subroutine adjoint_transform(K, Omega, dt)
        real(DP), intent(inout) :: K(DIM, DIM)
        real(DP), intent(in) :: Omega(DIM, DIM), dt
        real(DP) :: H(DIM, DIM), H_inv(DIM, DIM), tmp(DIM, DIM)
        call matrix_exp(Omega * dt, H)
        call invert_matrix(H, H_inv)
        tmp = matmul(H, K)
        K = matmul(tmp, H_inv)
    end subroutine adjoint_transform
end module pkgf_math_utils
