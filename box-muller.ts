// Box–Muller implementation
// https://en.wikipedia.org/wiki/Box%e2%80%93Muller_transform

const PRECISION = 1e9; // TODO: can remove?
const PI2 = Math.PI * 2;

/**
 * GenerateGaussianNoise
 * @param {number} mu
 * @param {number} sigma
 * @returns {number}
 */
export function generateGaussianNoise(mu: number, sigma: number): number {
  const u1 = Math.random();
  const u2 = Math.random();

  const z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(PI2 * u2);
  const z1 = Math.sqrt(-2.0 * Math.log(u1)) * Math.sin(PI2 * u2); // TODO: can remove?

  return z0 * sigma + mu;
}
