"""Exact certificates for the physical track witness and branch-moment controls."""
from fractions import Fraction as Q


def bump_bounds():
    # B=(1-q^2)^4, q=2z-1. (B')^2=256*t(1-t)^6, t=q^2.
    # Max at t=1/7. Negative curvature magnitude <= 32
    # for t<=1/7; positive magnitude max at 3/7.
    assert 256*Q(1,7)*Q(6,7)**6 < 16
    assert Q(4,7)**2 * 2 < 1
    # B'' at midpoint is exactly -32; joins at the endpoints are C^3.
    print("bump bounds: B<=1, |B'|<=4, |B''|<=32")


def track():
    e=Q(1,200000)
    du,dt=Q(1,16),Q(5,256)
    u_min,u_max=Q(3,2),Q(25,16)
    x_min,x_max=Q(5,4),Q(11,8)
    assert u_min*u_min/2>1 and u_max*u_max/2<x_min
    assert (x_min-1)**2/4==Q(1,64)
    assert (x_max-1)**2/4==Q(1,64)+dt
    G,P=Q(3,8),Q(3)
    assert (x_max-1)**2 * x_max/2<G*G
    assert (x_max-1)/2+x_max<P
    k1=4/du+4*G/dt
    k2=32/du**2+32*G/(du*dt)+32*G**2/dt**2+4*P/dt
    assert k1<141 and k2<31000
    delta=e*141
    L=1/(1-delta)
    dH=3*(L+1)*delta*L+G*L**3*31000*e
    assert dH<Q(1,8)
    # A strict source-strip geometry guarantee: |S_v| and its transformed
    # value both exceed 1/100 throughout support.
    sv_sq=(x_min-1)**2/4*(2*x_min-u_max*u_max)
    assert sv_sq>Q(1,100)**2*(1+delta)**2
    print("global positivity: Hessian perturbation < 1/8; geometric jets match for |c|<1/100")

    # Midpoint x=1+sqrt(13/128); both normalized bumps equal one.
    assert (Q(1,64)+dt/2)*4==Q(13,128)
    xlo,xhi=Q(131,100),Q(33,25)
    assert (xlo-1)**2<Q(13,128)<(xhi-1)**2
    u=Q(49,32)
    assert u==(u_min+u_max)/2
    assert 2*xlo-u*u>Q(1,4) and 2*xhi-u*u<1
    sv_lo=(xlo-1)/4
    grad_sq_lo=(xlo-1)**2*xlo/2
    lap_K_mag_lo=32*(1/du**2+grad_sq_lo/dt**2)
    same_point_gap=e*sv_lo*lap_K_mag_lo/2-e-e**2/2
    assert same_point_gap>Q(1,400)
    assert Q(3,2)-xlo==Q(19,100)
    assert same_point_gap>Q(1,100)*Q(19,100)
    assert xhi+e+e**2/2<Q(7,5)
    print("same-point convergence deficit >1/400 and >1% at a certified physical lens point")

    # The translated interval on the chosen fibre has same length 2v,
    # but an increased first moment 2*v*e>e since v>1/2.
    assert 2*xlo-u*u>Q(1,4)
    print("fibre length unchanged, first transverse moment strictly changes at midpoint")


def radial():
    a,b=Q(13,16),Q(19,16)
    t=(a-1)**2/4
    assert t==(b-1)**2/4==Q(9,1024)
    z=(t-Q(1,256))/(Q(1,64)-Q(1,256))
    assert z==Q(5,12)
    shift=Q(1,200000)*256*z**4*(1-z)**4
    A0=b-a
    A1=(b*b-a*a)/2
    c,d=a+shift,b+shift
    C0=d-c
    C1=(d*d-c*c)/2
    assert C0==A0 and C1-A1==shift*A0>0
    assert (C1/C0-C0/2,C1/C0+C0/2)==(c,d)
    # Angular derivative remainder: mean cos^4/3 = 1/8.
    assert Q(3,8)/3==Q(1,8)
    print("circular witness: aligned area exactly preserved; inverse endpoints recovered from first radial moment")


if __name__=="__main__":
    bump_bounds()
    track()
    radial()
    print("PASS")
