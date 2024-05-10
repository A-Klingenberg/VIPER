import logging
from abc import ABCMeta, abstractmethod
from typing import List

import ConfigManager
from modules.wrappers.RosettaWrapper import REBprocessor

cm = ConfigManager.ConfigManager.get_instance


class ResSelectionStrategy(metaclass=ABCMeta):
    """
    Defines behavior a selection strategy needs to implement. Namely, given a list of nodes, select a subset as a
    peptide candidate.
    """

    def __init__(self):
        self.max_length = cm().get("peptide_generator.max_length")
        self.linker = cm().get("peptide_generator.linker")
        self.linking_force_length_limit = cm().get("peptide_generator.linking_force_length_limit")
        self.length_damping = cm().get("peptide_generator.length_damping")
        self.ld_min_length = abs(cm().get("peptide_generator.length_damping_min_length"))
        self.ld_max_length = int(abs(cm().get("peptide_generator.length_damping_max_length")))
        if self.ld_min_length >= self.ld_max_length:
            logging.warning(f"ld_min ({self.ld_min_length}) is bigger than ld_max ({self.ld_max_length})!")
            if not cm().get("permissive"):
                raise ValueError(f"ld_min ({self.ld_min_length}) is bigger than ld_max ({self.ld_max_length})!")
        self.ld_initial_mult = abs(cm().get("peptide_generator.length_damping_initial_mult"))
        self.ld_final_mult = abs(cm().get("peptide_generator.length_damping_final_mult"))
        self.ld_linear_stepping = cm().get("peptide_generator.length_damping_linear_stepping", None)
        if m := cm().get("peptide_generator.length_damping_mode"):
            if self.ld_initial_mult > self.ld_final_mult:
                if m.upper() == "QUADRATIC":
                    self.damping_func = self._damping_factor_quadratic_penalty
                elif m.upper() == "LINEAR":
                    self.damping_func = self._damping_factor_linear_penalty
                else:
                    logging.warning(f"Could not parse damping mode {m} for SelectionStrategy. "
                                    f"Using default (QUADRATIC)...")
                    self.damping_func = self._damping_factor_quadratic_penalty
            else:
                if m.upper() == "QUADRATIC":
                    self.damping_func = self._damping_factor_quadratic_bonus
                elif m.upper() == "LINEAR":
                    self.damping_func = self._damping_factor_linear_bonus
                else:
                    logging.warning(f"Could not parse damping mode {m} for SelectionStrategy. "
                                    f"Using default (QUADRATIC)...")
                    self.damping_func = self._damping_factor_quadratic_bonus
        if self.ld_linear_stepping is not None:
            if self.damping_func(self.ld_max_length) != self.ld_final_mult:
                logging.warning(
                    f"Your configured settings for the linear damping factor yield a function that is not "
                    f"continuous. There will be a jump when reaching the length_damping_max_length!")
        else:
            self.ld_linear_stepping = (self.ld_initial_mult - self.ld_final_mult) / (
                    self.ld_max_length - self.ld_min_length)

    @abstractmethod
    def reduce(self, from_chain: str, to_chain: str, nodes: List[REBprocessor.Node]) -> List[REBprocessor.Node]:
        raise NotImplementedError("Trying to use ResSelectionStrategy:reduce() from abstract base class!")

    def _damping_factor_quadratic_penalty(self, peptide_length: int) -> float:
        """
        Returns the quadratic damping factor.
        The damping factor is the value of a parabola open to the bottom, centered on the minimum length
        with a maximum of the initial bonus. Any length before the minimum length is equal to the initial
        multiplier. Any length after and including the max length is equal to the final multiplier.

                        ∧ (damping factor)
                        |
        initial bonus   |---*______
                        |          ‾‾‾‾----___
                        |                     ‾‾--_
        max penalty     |                          ‾*---------
                        |
        ----------------+---|-----------------------|----------> (peptide length)
                        |   min length               max length

        :param peptide_length: Length of the peptide that the damping factor should be calculated for
        :return: Quadratic damping factor
        """
        return max(self.ld_final_mult,
                   (-1 * (max(peptide_length, self.ld_min_length) - self.ld_min_length) ** 2) *
                   max(0, (self.ld_initial_mult - self.ld_final_mult) /
                       (max(self.ld_min_length + 0.1, self.ld_max_length) - self.ld_min_length) ** 2) +
                   self.ld_initial_mult)

    def _damping_factor_quadratic_bonus(self, peptide_length: int) -> float:
        """
        Returns the quadratic damping factor.
        The damping factor is the value of a parabola centered on the minimum length with a maximum of the initial
        bonus. Any length before the minimum length is equal to the initial multiplier. Any length after and
        including the max length is equal to the final multiplier.

                        ∧ (damping factor)
                        |
        initial bonus   |                          _*---------
                        |                     __--‾
                        |          ____----‾‾‾
        max penalty     |---*‾‾‾‾‾‾
                        |
        ----------------+---|-----------------------|----------> (peptide length)
                        |   min length               max length

        :param peptide_length: Length of the peptide that the damping factor should be calculated for
        :return: Quadratic damping factor
        """
        return max(self.ld_initial_mult,
                   (min(self.ld_final_mult, ((max(peptide_length, self.ld_min_length) - self.ld_min_length) ** 2) *
                        max(0.0, (self.ld_final_mult - self.ld_initial_mult) /
                            ((max(self.ld_min_length + 0.1, self.ld_max_length) - self.ld_min_length) ** 2)) +
                        self.ld_initial_mult)))

    def _damping_factor_linear_penalty(self, peptide_length: int) -> float:
        """
        Returns the linear damping factor.
        The damping factor is the initial bonus, discounted by the stepping for every residue above the min length,
        up to the max length.

                        ∧ (damping factor)
                        |
        initial bonus  -|---*___                (The slope equals the stepping)
                        |       ‾‾‾---___
                        |                ‾‾‾---___
        max penalty     |                         ‾‾‾*---------
                        |
        ----------------+---|------------------------|----------> (peptide length)
                        |   min length               max length


        :param peptide_length: Length of the peptide that the damping factor should be calculated for
        :return: Linear damping factor
        """
        return max(max(0, self.ld_final_mult),
                   self.ld_initial_mult - min(max(0, peptide_length - self.ld_min_length),
                                              self.ld_max_length) * self.ld_linear_stepping)

    def _damping_factor_linear_bonus(self, peptide_length: int) -> float:
        """
        Returns the linear damping factor.
        The damping factor is the initial bonus, discounted by the stepping for every residue above the min length,
        up to the max length.

                        ∧ (damping factor)
                        |
        initial bonus   |                         ___*--------
                        |                ___---‾‾‾
                        |       ___---‾‾‾
        max penalty    -|---*‾‾‾                (The slope equals the stepping)
                        |
        ----------------+---|------------------------|----------> (peptide length)
                        |   min length               max length


        :param peptide_length: Length of the peptide that the damping factor should be calculated for
        :return: Linear damping factor
        """
        return min(max(0, self.ld_final_mult),
                   self.ld_initial_mult - min(max(0, peptide_length - self.ld_min_length),
                                              self.ld_max_length) * self.ld_linear_stepping)
