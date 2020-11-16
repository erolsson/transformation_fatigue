import numpy as np

from transformation_fatigue.materials.materials import SS2506


def m_stress(a1, a2, a3, Mss):
    p = np.array([-65, -230, -115.7])
    s_vm = np.array([736, 675, 659])
    r = np.array([580, 677, 589])
    return -3*a1*p + a2*s_vm**2/3 + 2./27*a3*r**3 + Mss


def main():
    # smooth center, smooth edge, notched center
    print(m_stress(0.02, 1e-4, 0, -56.05272438731849))
    print(m_stress(0.0, 1e-4, 6e-7, -56.05272438731849))
    print(m_stress(0.03, 2e-4, 0, -75))
    c = 0.007989
    k = SS2506.k_1 + SS2506.k_2*c
    Ms = SS2506.Ms_1 + SS2506.Ms_2*c
    print(-np.log(0.2105)/k - Ms + 22.)


if __name__ == '__main__':
    main()
