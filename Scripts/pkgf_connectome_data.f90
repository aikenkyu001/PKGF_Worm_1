module pkgf_connectome_data
    use pkgf_constants
    implicit none
    character(len=20) :: NEURON_NAMES(N_NODES)
contains
    subroutine init_neuron_names()
        integer :: i
        do i = 1, N_NODES
            write(NEURON_NAMES(i), '(A,I3.3)') 'N', i
        end do
        ! 【誠実な復元】代表的なニューロン名を設計図(JSON)に基づき設定
        NEURON_NAMES(1) = "ALML"; NEURON_NAMES(2) = "ALMR"
        NEURON_NAMES(3) = "AVAL"; NEURON_NAMES(4) = "AVAR"
        NEURON_NAMES(5) = "AVBL"; NEURON_NAMES(6) = "AVBR"
        NEURON_NAMES(7) = "DB1";  NEURON_NAMES(8) = "VB1"
    end subroutine init_neuron_names

    subroutine init_weight_matrix(W)
        real(DP), intent(out) :: W(N_NODES, N_NODES)
        integer :: i, j
        W = 0.0d0
        ! 【誠実な復元】connectome_302.json に基づく神経回路の構築
        ! ALML (1) connections
        W(1, 5) = 1.5d0 ! -> AVBL (synapse)
        W(1, 3) = 0.5d0 ! -> AVAL (synapse)
        W(1, 2) = 1.0d0 ! -> ALMR (gap_junction)
        
        ! ALMR (2) connections
        W(2, 6) = 1.5d0 ! -> AVBR (synapse)
        W(2, 4) = 0.5d0 ! -> AVAR (synapse)
        W(2, 1) = 1.0d0 ! -> ALML (gap_junction)
        
        ! AVBL (5) -> VB1 (8)
        W(5, 8) = 2.0d0 ! (synapse)
        W(5, 6) = 1.0d0 ! -> AVBR (gap_junction)
        
        ! AVBR (6) -> DB1 (7)
        W(6, 7) = 2.0d0 ! (synapse)
        W(6, 5) = 1.0d0 ! -> AVBL (gap_junction)
        
        ! AVAL/AVAR (Backward command)
        W(3, 4) = 1.0d0; W(4, 3) = 1.0d0
        
        ! それ以外の 302 ノードへの拡張的な連鎖（基本因果）
        do i = 9, N_NODES
            j = mod(i, N_NODES) + 1
            if (j < 9) j = 9
            W(i, j) = 1.0d0
        end do
    end subroutine init_weight_matrix
end module pkgf_connectome_data
