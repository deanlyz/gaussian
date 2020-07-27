import { generateGaussianNoise } from "./box-muller.ts";

export class Gaussian {
  readonly mean: number;
  readonly variance: number;
  readonly standardDeviation: number;

  constructor(mean: number, variance: number) {
    if (variance <= 0) {
      throw new Error(`Variance must be > 0 (but was: ${variance})`);
    }

    this.mean = mean;
    this.variance = variance;
    this.standardDeviation = Math.sqrt(variance);
  }

  private static fromPrecisionMean(precision: number, precisionMean: number) {
    return new Gaussian(precisionMean / precision, 1 / precision);
  }

  /**
     * Probability density function
     * @param {number} x
     * @return {number}
     */
  pdf(x: number) {
    const m = this.standardDeviation * Math.sqrt(2 * Math.PI);
    const e = Math.exp(-Math.pow(x - this.mean, 2) / (2 * this.variance));
    return e / m;
  }

  /**
     * Cumulative density function
     * @param {number} x
     * @return {number}
     */
  cdf(x: number): number {
    return 0.5 *
      erfc(-(x - this.mean) / (this.standardDeviation * Math.sqrt(2)));
  }

  /**
     * Percent point function
     * @param {number} x
     * @return {number}
     */
  ppf(x: number): number {
    return this.mean - this.standardDeviation * Math.sqrt(2) * ierfc(2 * x);
  }

  /**
     * Product distribution of this and d (scale for constant)
     * @param {number} d
     * @return {Gaussian}
     */
  mul(d: number | Gaussian): Gaussian {
    if (typeof d === "number") {
      return this.scale(d);
    }

    const precision = 1 / this.variance;
    const dPrecision = 1 / d.variance;
    return Gaussian.fromPrecisionMean(
      precision + dPrecision,
      precision * this.mean + dPrecision * d.mean,
    );
  }

  /**
     * Quotient distribution of this and d (scale for constant)
     * @param {number | Gaussian} d
     * @return {Gaussian}
     */
  div(d: number | Gaussian): Gaussian {
    if (typeof d === "number") {
      return this.scale(1 / d);
    }

    const precision = 1 / this.variance;
    const dPrecision = 1 / d.variance;
    return Gaussian.fromPrecisionMean(
      precision - dPrecision,
      precision * this.mean - dPrecision * d.mean,
    );
  }

  /**
     * Addition of this and d
     * @param {number | Gaussian} d
     * @return {Gaussian}
     */
  add(d: Gaussian): Gaussian {
    return new Gaussian(this.mean + d.mean, this.variance + d.variance);
  }

  /**
     * Subtraction of this and d
     * @param {number | Gaussian} d
     * @return {Gaussian}
     */
  sub(d: Gaussian): Gaussian {
    return new Gaussian(this.mean - d.mean, this.variance + d.variance);
  }

  /**
     * Scale this by constant c
     * @param {number} c
     * @return {Gaussian}
     */
  scale(c: number): Gaussian {
    return new Gaussian(this.mean * c, this.variance * c * c);
  }

  /**
     * Generate [num] random samples
     * @param {number} num
     * @return {number[]}
     */
  random(num: number): number[] {
    const mean = this.mean;
    const std = this.standardDeviation;
    return Array(num).fill(0).map(() => generateGaussianNoise(mean, std));
  }
}

/**
 * Create a Gaussian
 * @param {number} mean
 * @param {number} variance
 * @return {Gaussian}
 */
export function gaussian(mean: number, variance: number): Gaussian {
  return new Gaussian(mean, variance);
}

/**
 * Complementary error function
 * From Numerical Recipes in C 2e p221
 * @param {number} x
 * @return {number}
 */
export function erfc(x: number): number {
  const z = Math.abs(x);
  const t = 1 / (1 + z / 2);
  const r = t *
    Math.exp(
      -z * z - 1.26551223 +
        t *
          (1.00002368 +
            t *
              (0.37409196 +
                t *
                  (0.09678418 +
                    t *
                      (-0.18628806 +
                        t *
                          (0.27886807 +
                            t *
                              (-1.13520398 +
                                t *
                                  (1.48851587 +
                                    t * (-0.82215223 + t * 0.17087277)))))))),
    );

  return x >= 0 ? r : 2 - r;
}

/**
 * Inverse complementary error function
 * From Numerical Recipes 3e p265
 * @param {number} x
 * @return {number}
 */
export function ierfc(x: number): number {
  if (x >= 2) return -100;
  if (x <= 0) return 100;

  const xx = x < 1 ? x : 2 - x;
  const t = Math.sqrt(-2 * Math.log(xx / 2));

  let r = -0.70711 *
    ((2.30753 + t * 0.27061) / (1 + t * (0.99229 + t * 0.04481)) - t);

  for (let i = 0; i < 2; i++) {
    const e = erfc(r) - xx;
    r += e / (1.12837916709551257 * Math.exp(-(r * r)) - r * e);
  }

  return x < 1 ? r : -r;
}
