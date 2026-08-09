#!/usr/bin/env python3

import unittest

from scripts.thread_broker_oracle import oracle_winner_is_bracketed


class OracleBoundaryTests(unittest.TestCase):
    def test_interior_winner_is_bracketed(self):
        self.assertTrue(
            oracle_winner_is_bracketed(
                {"pin-1": 9.0, "pin-2": 8.0, "pin-4": 10.0}, "pin-2"
            )
        )

    def test_low_boundary_winner_is_not_an_oracle(self):
        self.assertFalse(
            oracle_winner_is_bracketed(
                {"pin-2": 8.0, "pin-4": 10.0, "pin-8": 12.0}, "pin-2"
            )
        )

    def test_high_boundary_winner_is_not_an_oracle(self):
        self.assertFalse(
            oracle_winner_is_bracketed(
                {"pin-2": 12.0, "pin-4": 10.0, "pin-8": 8.0}, "pin-8"
            )
        )


if __name__ == "__main__":
    unittest.main()
