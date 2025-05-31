"""
1D Line Plot
============
"""

import numpy as np
import matplotlib.pyplot as plt

from reservoirflow.plots.plot import Plot


class Plot1D(Plot):
    """1D Line Plot class for visualizing 1D data.

    Parameters
    ----------
    Plot : _type_
        _description_
    """

    name = "1D Line Plot"

    def __init__(
        self,
        nrows=1,
        ncols=1,
        verbose=False,
    ):
        """Create 1D Line Plot object.

        Parameters
        ----------
        nrows : int
            number of rows in the plot.
        ncols : int
            number of columns in the plot.
        verbose : bool, optional
            print verbose output, by default False
        """
        super().__init__(verbose)
        self.nrows = nrows
        self.ncols = ncols

    def __set_axis_labels(self, axs):
        for i, ax in enumerate(axs.ravel()):
            if i >= 6:
                ax.set_xlabel("x")
            if i in [0, 3, 6]:
                ax.set_ylabel("p")

    def __add_legends(self, fig, N):
        # https://stackoverflow.com/a/59393045/11549398
        labels_handles = {
            label: handle
            for ax in fig.axes
            for handle, label in zip(*ax.get_legend_handles_labels())
        }

        fig.legend(
            labels_handles.values(),
            labels_handles.keys(),
            loc="upper center",
            bbox_to_anchor=(0.5, 0.03),
            ncol=6,  # ncol=N,
        )

    def __get_lims_ticks(self, x, Y):

        xmin = x.min()
        xmax = x.max()
        xmax_ = round(xmax, 1)
        xstep = (xmax_ - xmin) / 4
        xmin_ = xmin - xstep * 0.25
        xlim = (xmin_, xmax * 1.1)
        xticks = np.linspace(xmin, xmax_, 5)

        ymin = Y.min()
        ymax = Y.max()
        n = 0 if ymax < 2 else -3
        ymax_ = round(ymax, n)
        ystep = (ymax_ - ymin) / 4
        ymin_ = ymin - ystep * 0.25
        ylim = (ymin_, ymax_)
        if ymax_ > ymax:
            ystep = ystep if ymax < 2 else 0
        m = 5 if ystep == 0 else 6
        yticks = np.linspace(ymin, ymax_ + ystep, m)

        return (xlim, xticks), (ylim, yticks)

    def plot(
        self,
        id=None,
    ):

        def get_sampled_array_combined_distribution(input_array, num_samples):
            """
            Selects an exact number of values from a 1D NumPy array using a
            combination of logarithmic (geomspace-like) and linear distributions.

            - The first segment of selected samples approximates a logarithmic distribution.
            - The second segment of selected samples approximates a linear distribution.

            Guarantees:
            - The exact number of samples specified by `num_samples` is returned.
            - Index 0 and the last index (len(input_array)-1) are always included,
            unless `num_samples` is too small to cover all unique indices (e.g.,
            if num_samples is 1 and the array has many elements, only the first is returned).
            - The desired distribution type (log then linear) is approximated.

            Parameters:
            -----------
            input_array : np.ndarray
                The 1D NumPy array from which to select values.
            num_samples : int
                The exact desired number of values to select. Must be >= 1.

            Returns:
            --------
            np.ndarray
                A 1D array containing exactly `num_samples` selected values.
            """
            import math

            if input_array.ndim != 1:
                raise ValueError("Input array must be 1-dimensional.")

            num_elements = len(input_array)

            if num_elements == 0:
                return np.array([])  # Handle empty input array

            if num_samples < 1:
                raise ValueError("num_samples must be at least 1.")

            max_index = num_elements - 1

            # --- Handle base cases for very small num_samples ---
            if num_samples == 1:
                return np.array([input_array[0]])
            if num_samples == 2:
                # Guarantee first and last elements for num_samples = 2
                return np.array([input_array[0], input_array[-1]])

            # --- Handle case where requested samples exceed available unique elements ---
            if num_samples > num_elements:
                # Generate all unique indices available in the array
                all_unique_indices_available = np.arange(num_elements)
                # Pad the result by repeating the last element to reach num_samples
                padded_indices = np.concatenate(
                    (
                        all_unique_indices_available,
                        np.full(
                            num_samples - num_elements, all_unique_indices_available[-1]
                        ),
                    )
                )
                return input_array[padded_indices]

            # --- Generate a sufficiently large pool of unique candidate indices ---
            # We generate more candidates than strictly needed to ensure we have enough
            # unique integer points after rounding and deduplication, covering the desired distribution.
            # The multiplier (e.g., 5) can be adjusted for denser initial candidate pools.
            candidate_pool_size = max(
                num_samples * 5, 200
            )  # Ensure a minimum pool size for better distribution approximation

            # Determine the number of raw points to generate for each segment in the candidate pool
            num_geom_candidates = math.ceil(candidate_pool_size * 2 / 3)
            num_lin_candidates = candidate_pool_size - num_geom_candidates

            # Ensure non-zero candidate counts for distribution functions to work
            if num_geom_candidates == 0:
                num_geom_candidates = 1
            if num_lin_candidates == 0:
                num_lin_candidates = 1

            # Define the transition point in the array's index range (approx. 2/3 of max_index)
            split_index_value = max_index * (2 / 3.0)

            # Generate floating-point indices for the geomspace-like part (from ~0 to split_index_value)
            candidate_geom_indices_float = np.array([])
            if num_geom_candidates > 0:
                # np.geomspace needs positive start/end. We use a small offset (1e-9) to simulate 0.
                # The range is [offset, split_index_value + offset]
                geom_start_val = 1e-9
                geom_end_val = max(
                    1e-9, split_index_value
                )  # Ensure positive end for geomspace
                candidate_geom_indices_float = (
                    np.geomspace(
                        geom_start_val,
                        geom_end_val + geom_start_val,
                        num_geom_candidates,
                    )
                    - geom_start_val
                )
                # Clip to ensure bounds [0, split_index_value]
                candidate_geom_indices_float = np.clip(
                    candidate_geom_indices_float, 0, split_index_value
                )

            # Generate floating-point indices for the linspace part (from split_index_value to max_index)
            candidate_lin_indices_float = np.array([])
            if num_lin_candidates > 0:
                candidate_lin_indices_float = np.linspace(
                    split_index_value, max_index, num_lin_candidates
                )
                # Clip to ensure bounds [split_index_value, max_index]
                candidate_lin_indices_float = np.clip(
                    candidate_lin_indices_float, split_index_value, max_index
                )

            # Combine all candidate indices, add explicit 0 and max_index to guarantee their presence
            all_candidate_indices_float = np.concatenate(
                (
                    [0.0],  # Ensure 0 is always included
                    candidate_geom_indices_float,
                    candidate_lin_indices_float,
                    [float(max_index)],  # Ensure max_index is always included
                )
            )

            # Round to integer, then get a unique and sorted list of candidate indices from the large pool
            unique_sorted_indices_pool = np.unique(
                np.round(all_candidate_indices_float)
            ).astype(int)

            # Ensure all indices are strictly within valid array bounds
            unique_sorted_indices_pool = np.clip(
                unique_sorted_indices_pool, 0, max_index
            )

            # --- Final selection to guarantee EXACT `num_samples` indices ---
            # We select `num_samples` points by evenly spacing through the `unique_sorted_indices_pool`.
            # This method effectively picks `num_samples` points that best represent the overall distribution
            # of the larger candidate pool.

            if len(unique_sorted_indices_pool) < num_samples:
                # This case happens if max_index is very small, and even the large pool can't yield
                # enough unique integer indices. In this scenario, we must repeat points.
                # We fill by repeating the last unique index we could get.
                final_selected_indices = np.concatenate(
                    (
                        unique_sorted_indices_pool,
                        np.full(
                            num_samples - len(unique_sorted_indices_pool),
                            unique_sorted_indices_pool[-1],
                        ),
                    )
                )
            else:
                # Calculate indices to select from the unique_sorted_indices_pool
                # np.linspace here creates `num_samples` evenly spaced integer positions within the pool's ranks.
                selection_indices_in_pool = np.linspace(
                    0, len(unique_sorted_indices_pool) - 1, num_samples
                ).astype(int)
                final_selected_indices = unique_sorted_indices_pool[
                    selection_indices_in_pool
                ]

            # Return the values from the input array corresponding to the final selected indices, sorted.
            return np.sort(final_selected_indices)

        def get_sampled_array_flexible(input_array, num_samples):
            """
            Selects a flexible number of values from a 1D NumPy array.
            The selection ensures the first (index 0) and last (index len-1)
            elements are always included, with intermediate elements
            being approximately evenly spaced across the array.

            Parameters:
            -----------
            input_array : np.ndarray
                The 1D NumPy array from which to select values.
            num_samples : int
                The desired number of values to select.
                - Must be >= 1.
                - If 1, only the first element will be returned (repeated `num_samples` times if array has 1 element).
                - If 2 or more, both the first and last elements are guaranteed to be included.

            Returns:
            --------
            np.ndarray
                A 1D array containing the selected values.
                If `num_samples` is greater than the array's actual length,
                some values in the output array might be duplicates.
            """

            if input_array.ndim != 1:
                raise ValueError("Input array must be 1-dimensional.")

            num_elements = len(input_array)

            if num_elements == 0:
                return np.array([])  # Handle empty input array

            if num_samples < 1:
                raise ValueError(
                    "The number of samples (num_samples) must be at least 1."
                )

            if num_elements == 1:
                # If the input array has only one element, repeat it num_samples times
                return np.array([input_array[0]] * num_samples)

            # Calculate the maximum possible index for the current array length
            max_index = num_elements - 1

            # Use np.linspace to generate evenly spaced floating-point indices
            # between 0 (first element) and max_index (last element).
            # np.linspace handles cases where num_samples is 1 or 2 correctly.
            selected_indices_float = np.linspace(0, max_index, num_samples)

            # Round the floating-point indices to the nearest integer
            # and convert them to integer type for array indexing.
            selected_indices = np.round(selected_indices_float).astype(int)

            # Use np.clip as a safety measure to ensure indices are within valid bounds.
            # This is typically not needed with np.linspace and np.round for this range,
            # but it adds robustness.
            selected_indices = np.clip(selected_indices, 0, max_index)

            # Select the values from the input array using the calculated indices
            return selected_indices

        fig, axs = plt.subplots(
            self.nrows,
            self.ncols,
            figsize=(10, 6),
            sharey=True,
            sharex=True,
        )

        plt.subplots_adjust(hspace=0.3, wspace=0.2)
        alpha = 1
        # tsteps = [0, 1, 5, 10, 40, 50, 60, 90, 100]

        # assert len(tsteps) == self.nrows * self.ncols, "tsteps are not compatible"

        self.__set_axis_labels(axs)

        if id is None:
            names = list(self.Data.keys())
            N = len(names)
        else:
            if isinstance(id, list):
                lst = list(self.Data.keys())
                names = [lst[i] for i in id]
                N = len(names)
            else:
                names = [list(self.Data.keys())[id]]
                N = 1

        for n in range(N):
            name = names[n]
            lstyle = "-" if "ana" in name.lower() else "--"
            X = self.Data[name][0]
            Y = self.Data[name][1]
            t = X[:, 0, 0]
            x = X[0, :, 1]

            tsteps = get_sampled_array_combined_distribution(t, self.nrows * self.ncols)

            assert len(tsteps) == self.nrows * self.ncols, "tsteps are not compatible"

            for i, ax in enumerate(axs.ravel()):
                tstep = tsteps[i]
                ax.plot(x, Y[tstep, :], label=name, linestyle=lstyle, alpha=alpha)
                ax.tick_params(top=True, right=True)
                ax.set_title(f"t={t[tstep]}")
                ax.grid(True)
                if i == 0 and n == 0:
                    (xlim, xticks), (ylim, yticks) = self.__get_lims_ticks(x, Y)
                    ax.set_xlim(xlim)
                    ax.set_ylim(ylim)
                    ax.set_xticks(xticks)
                    ax.set_yticks(yticks)

        self.__add_legends(fig, N)
        return plt.show()

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":
    plot = Plot1D()
    print(plot)
