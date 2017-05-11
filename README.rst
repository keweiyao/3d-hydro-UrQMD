hic-ebe-3d
==================

vhlle freezeout data:

Old data:

.. math::
  \tau, x, y, \eta, d\sigma_0, d\sigma_1, d\sigma_2, d\sigma_3, u^0, u^1, u^2, u^3, T_c, \mu_b, \mu_q, \mu_s, \pi_{00}, \pi_{01}, \pi_{02}, \pi_{03}, \pi_{11}, \pi_{12}, \pi_{13}, \pi_{22}, \pi_{23}, \pi_{33}, \pi_{33}, \Pi, dV_{eff}=d\sigma \cdot u


.. csv-table:: New data (everything in cartesian coordinates)
   :header: "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11-15"
   :widths: 3,3,3, 3,3,3, 3,3,3, 3,3,3

   :math:`x^0`, :math:`x^1`, :math:`x^2`, :math:`x^3`, :math:`d\sigma_0`, :math:`d\sigma_1`, :math:`d\sigma_2`, :math:`d\sigma_3`, :math:`v_1`, :math:`v_2`, :math:`v_3`, ":math:`\pi_{ij = (11, 12, 13, 22, 23)}`"


Note that sampler needs:

.. math::
   x = [x^0, x^1, x^2, x^3] & &\\
   d\sigma^{\mu} = [d\sigma_0, -d\sigma_1, -d\sigma_2, -d\sigma_3] & &\\
   v = [v_1, v_2, v_3] & &\\
   \pi_{ij} \textrm{ or } \pi^{ij} \textrm{ dictionary } & &\\

and make sure that :code:`Sampler.volume` is positive.


