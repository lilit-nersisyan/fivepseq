import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


def plot_frames_term(fivepseqcounts, region, pdf_file):
    span_size = fivepseqcounts.span_size

    with PdfPages(pdf_file) as pdf:
        fig = plt.figure()
        plt.title(("Peak at: " + str(np.argmax(fivepseqcounts.get_meta_count_series(region)) - span_size)))
        plt.plot(fivepseqcounts.get_meta_count_series(region))
        plt.xlabel("Position at %s" % region)
        plt.close()
        pdf.savefig(fig)

        fig = plt.figure()
        plt.title(("Peak at: " + str(np.argmax(fivepseqcounts.get_meta_count_series(region)) - span_size)))
        plt.xlabel("Position at %s" % region)
        counts = fivepseqcounts.get_meta_count_series(region)
        frames = fivepseqcounts.extract_frame_count_vectors(counts)
        plt.plot(np.arange(1, len(counts) + 1), counts, '0.75',
                 np.arange(1, len(counts) + 1, 3), frames[0], 'ro',
                 np.arange(2, len(counts) + 1, 3), frames[1], 'go',
                 np.arange(3, len(counts) + 1, 3), frames[2], 'bo')
        plt.close()
        pdf.savefig(fig)

        fig = plt.figure()
        frames = fivepseqcounts.extract_frame_count_vectors(fivepseqcounts.get_meta_count_series(region))

        plt.plot(frames[0], 'r', frames[1], 'g', frames[2], 'b')
        plt.xlabel("Frame position at %s" % region)
        plt.ylabel("Number of 5' mappings")

        plt.close()
        pdf.savefig(fig)
