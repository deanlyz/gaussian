// Tests based on values from Wolfram Alpha.

import { assert, assertEquals, assertThrows } from "https://deno.land/std/testing/asserts.ts";

import { Gaussian } from "./gaussian.ts";

function assertDiffEquals(
    actual: number,
    expected: number,
    precision: number = 5,
): void {
  const actualDiff = Math.abs(actual - expected)
  const expectedDiff = Math.pow(10, -precision)
  return assert(actualDiff < expectedDiff) // TODO: add msg
}

Deno.test("test properties", () => {
  const d1 = new Gaussian(0, 1);
  assertEquals(d1.mean, 0);
  assertEquals(d1.variance, 1);
  assertEquals(d1.standardDeviation, 1);

  const d2 = new Gaussian(1, 4);
  assertEquals(d2.mean, 1);
  assertEquals(d2.variance, 4);
  assertEquals(d2.standardDeviation, 2);
});

Deno.test("test pdf", () => {
  const d = new Gaussian(0, 1);
  assertDiffEquals(d.pdf(-2), 0.053991)
  assertDiffEquals(d.pdf(-1), 0.241971)
  assertDiffEquals(d.pdf(0), 0.398942)
  assertDiffEquals(d.pdf(1), 0.241971)
  assertDiffEquals(d.pdf(2), 0.053991)
});

Deno.test('test cdf', () => {
  const d = new Gaussian(0, 1)
  assertDiffEquals(d.cdf(-1.28155), 0.1)
  assertDiffEquals(d.cdf(-0.67449), 0.25)
  assertDiffEquals(d.cdf(0), 0.5)
  assertDiffEquals(d.cdf(0.67449), 0.75)
  assertDiffEquals(d.cdf(1.28155), 0.9)
})

Deno.test('test ppf', () => {
  const d = new Gaussian(0, 1)
  assertDiffEquals(d.ppf(0.1), -1.28155)
  assertDiffEquals(d.ppf(0.25), -0.67449)
  assertDiffEquals(d.ppf(0.5), 0)
  assertDiffEquals(d.ppf(0.75), 0.67449)
  assertDiffEquals(d.ppf(0.9), 1.28155)
})

Deno.test('test mul', () => {
  const d = new Gaussian(0, 1).mul(new Gaussian(0, 1))
  assertEquals(d, new Gaussian(0, 0.5))
  assertEquals(new Gaussian(1, 1).mul(2), new Gaussian(2, 4))
})

Deno.test('test div', () => {
  const d = new Gaussian(1, 1).div(new Gaussian(1, 2))
  assertEquals(d, new Gaussian(1, 2))
  assertEquals(new Gaussian(1, 1).div(1/2), new Gaussian(2, 4))
})

Deno.test('test rejects non-positive variances', () => {
  assertThrows(() => new Gaussian(0, 0), Error, 'Variance must be > 0 (but was: 0)')
  assertThrows(() => new Gaussian(0, -1), Error, 'Variance must be > 0 (but was: -1)')
})

Deno.test('test add', () => {
  assertEquals(new Gaussian(1, 1).add(new Gaussian(1, 2)), new Gaussian(2, 3))
})

Deno.test('test sub', () => {
  assertEquals(new Gaussian(1, 1).sub(new Gaussian(1, 2)), new Gaussian(0, 3))
})

Deno.test('test scale', () => {
  assertEquals(new Gaussian(1, 1).scale(2), new Gaussian(2, 4))
})

Deno.test('test generate samples', () => {
  const outcomes = new Gaussian(0, 0.3).random(10)
  outcomes.forEach(o => assertEquals(typeof o, 'number'))
  assertEquals(outcomes.length, 10)
})

Deno.test('test generated sample distribution', () => {
  const size = 3e6
  const outcomes = new Gaussian(-1, 0.65).random(size)
  let mean = 0
  let variance = 0
  outcomes.forEach(n => mean += n)
  mean /= size

  outcomes.forEach(n => variance += Math.pow(n - mean, 2))
  variance /= size

  assertDiffEquals(mean, -1.0, 3)
  assertDiffEquals(variance, 0.65, 3)
})
