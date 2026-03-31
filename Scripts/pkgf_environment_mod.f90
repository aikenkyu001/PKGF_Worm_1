module pkgf_environment_mod
    use pkgf_constants
    implicit none
    type Food
        real(DP) :: pos(2), intensity
        logical :: consumed
    contains
        procedure :: consume => food_consume
    end type Food
    type Trace
        real(DP) :: pos(2), age
    end type Trace
    type EnvironmentalPKGF
        type(Food) :: foods(1)
        type(Trace) :: traces(200)
        integer :: n_traces
    contains
        procedure :: add_food => env_add_food
        procedure :: add_trace => env_add_trace
        procedure :: get_potential => env_get_potential
        procedure :: get_gradient => env_get_gradient
        procedure :: get_curvature => env_get_curvature
        procedure :: step => env_step
    end type EnvironmentalPKGF
contains
    function env_get_curvature(self, obs_p) result(F)
        class(EnvironmentalPKGF), intent(in) :: self
        real(DP), intent(in) :: obs_p(2)
        real(DP) :: F(2, 2), d_sq, diff(2), B_field
        F = 0.0d0
        B_field = 0.0d0
        if (.not. self%foods(1)%consumed) then
            diff = self%foods(1)%pos - obs_p
            d_sq = sum(diff**2) + ALPHA**2 ! 分母の安定化をシステム定数 ALPHA で実施
            B_field = self%foods(1)%intensity / d_sq
        end if
        F(1, 2) = B_field; F(2, 1) = -B_field
    end function env_get_curvature

    function food_consume(self, amount) result(actual)
        class(Food), intent(inout) :: self
        real(DP), intent(in) :: amount
        real(DP) :: actual
        if (self%consumed) then
            actual = 0.0d0
        else
            actual = amount ! 上限判定を廃止し、流動の必然に任せる
            self%intensity = self%intensity - actual
            if (self%intensity <= 0.0d0) then
                self%intensity = 0.0d0
                self%consumed = .true.
            end if
        end if
    end function food_consume

    subroutine env_add_food(self, pos, intensity)
        class(EnvironmentalPKGF), intent(inout) :: self
        real(DP), intent(in) :: pos(2), intensity
        self%foods(1)%pos = pos
        self%foods(1)%intensity = intensity
        self%foods(1)%consumed = .false.
    end subroutine env_add_food

    subroutine env_add_trace(self, pos)
        class(EnvironmentalPKGF), intent(inout) :: self
        real(DP), intent(in) :: pos(2)
        if (self%n_traces < 200) then
            self%n_traces = self%n_traces + 1
        else
            self%traces(1:199) = self%traces(2:200)
        end if
        self%traces(self%n_traces)%pos = pos
        self%traces(self%n_traces)%age = 1.0d0
    end subroutine env_add_trace

    function env_get_potential(self, obs_p) result(p)
        class(EnvironmentalPKGF), intent(in) :: self
        real(DP), intent(in) :: obs_p(2)
        real(DP) :: p, dist_sq
        integer :: i
        p = 0.0d0
        if (.not. self%foods(1)%consumed) then
            dist_sq = sum((self%foods(1)%pos - obs_p)**2) + ALPHA**2
            p = self%foods(1)%intensity / dist_sq
        end if
        do i = 1, self%n_traces
            dist_sq = sum((self%traces(i)%pos - obs_p)**2) + ALPHA**2
            p = p - (self%traces(i)%age / dist_sq)
        end do
    end function env_get_potential

    function env_get_gradient(self, obs_p) result(grad)
        class(EnvironmentalPKGF), intent(in) :: self
        real(DP), intent(in) :: obs_p(2)
        real(DP) :: grad(2), d_sq, diff(2)
        integer :: i
        grad = 0.0d0
        if (.not. self%foods(1)%consumed) then
            diff = self%foods(1)%pos - obs_p
            d_sq = sum(diff**2) + ALPHA**2
            grad = grad + (2.0d0 * self%foods(1)%intensity / (d_sq**2)) * diff
        end if
        do i = 1, self%n_traces
            diff = obs_p - self%traces(i)%pos
            d_sq = sum(diff**2) + ALPHA**2
            grad = grad + (2.0d0 * self%traces(i)%age / (d_sq**2)) * diff
        end do
    end function env_get_gradient

    subroutine env_step(self)
        class(EnvironmentalPKGF), intent(inout) :: self
        integer :: i
        ! 【浄化】根拠なき 0.98 倍の減衰を廃止。
        ! 本来は拡散方程式を解くべきだが、偽装するよりは減衰なし（保存）を選択する。
    end subroutine env_step
end module pkgf_environment_mod
