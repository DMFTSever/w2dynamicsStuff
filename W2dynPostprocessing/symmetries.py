import numpy as np


def time_rev_symm_LS_matsubara(arr):
    """Symmetrise array acording to time reversal symmetry:
       s1,s2 = up/down (+1/-1 for numerical purposes)
       A_{a,s1,b,s2}(iw)=s1*s2*A_{a,s1bar,b,s2bar}^*(-iw)
    """

    if arr.shape[-1] != 2:
        raise NotImplementedError("Only works for 2 spins for now")

    newarr = np.zeros_like(arr)
    newarr[..., :, 0, :, 0, :] = 0.5 * (arr[..., :, 0, :, 0, :] +
                                        arr[..., :, 1, :, 1, ::-1].conjugate())
    newarr[..., :, 0, :, 1, :] = 0.5 * (arr[..., :, 0, :, 1, :] -
                                        arr[..., :, 1, :, 0, ::-1].conjugate())
    newarr[..., :, 1, :, 0, :] = 0.5 * (arr[..., :, 1, :, 0, :] -
                                        arr[..., :, 0, :, 1, ::-1].conjugate())
    newarr[..., :, 1, :, 1, :] = 0.5 * (arr[..., :, 1, :, 1, :] +
                                        arr[..., :, 0, :, 0, ::-1].conjugate())

    return newarr[:]


def time_rev_symm_LS_tau(arr):
    """Symmetrise array acording to time reversal symmetry:
       s1,s2 = up/down (+1/-1 for numerical purposes)
       A_{a,s1,b,s2}(tau)=s1*s2*A_{a,s1bar,b,s2bar}^*(tau)
    """

    if arr.shape[-1] != 2:
        raise NotImplementedError("Only works for 2 spins for now")

    newarr = np.zeros_like(arr)
    newarr[..., :, 0, :, 0, :] = 0.5 * (arr[..., :, 0, :, 0, :] +
                                        arr[..., :, 1, :, 1, :].conjugate())
    newarr[..., :, 0, :, 1, :] = 0.5 * (arr[..., :, 0, :, 1, :] -
                                        arr[..., :, 1, :, 0, :].conjugate())
    newarr[..., :, 1, :, 0, :] = 0.5 * (arr[..., :, 1, :, 0, :] -
                                        arr[..., :, 0, :, 1, :].conjugate())
    newarr[..., :, 1, :, 1, :] = 0.5 * (arr[..., :, 1, :, 1, :] +
                                        arr[..., :, 0, :, 0, :].conjugate())

    return newarr[:]


def time_rev_symm_LS_scalar(arr):
    """Symmetrise array acording to time reversal symmetry:
       s1,s2 = up/down (+1/-1 for numerical purposes)
       A_{a,s1,b,s2}=s1*s2*A_{a,s1bar,b,s2bar}^*
    """

    if arr.shape[-1] != 2:
        raise NotImplementedError("Only works for 2 spins for now")

    newarr = np.zeros_like(arr)
    newarr[..., :, 0, :, 0] = 0.5 * (arr[..., :, 0, :, 0] +
                                     arr[..., :, 1, :, 1].conjugate())
    newarr[..., :, 0, :, 1] = 0.5 * (arr[..., :, 0, :, 1] -
                                     arr[..., :, 1, :, 0].conjugate())
    newarr[..., :, 1, :, 0] = 0.5 * (arr[..., :, 1, :, 0] -
                                     arr[..., :, 0, :, 1].conjugate())
    newarr[..., :, 1, :, 1] = 0.5 * (arr[..., :, 1, :, 1] +
                                     arr[..., :, 0, :, 0].conjugate())
