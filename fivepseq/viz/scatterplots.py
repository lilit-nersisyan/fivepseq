import matplotlib.pyplot as plt

def plot_frames(counts, pdf_file):
    plt.title(("Peak at: " + str(metagene[0].argmax() - (global_args.offset))))
    plt.grid(True, alpha=0.3)
    plt.step(np.linspace(-global_args.offset, global_args.offset, num=global_args.offset * 2),
             metagene[0], linewidth=0.5, label=bamname)
    plt.xticks(np.linspace(-global_args.offset, global_args.offset, num=global_args.offset * 2),
               labels, size="xx-small")
    plt.legend()
    plt.savefig(global_args.output_dir + "coverage_start/%s%s.pdf" % (bamname, norm))
    plt.close()

    plt.title(("Peak at: " + str(metagene[1].argmax() - (global_args.offset))))
    plt.grid(True, alpha=0.3)
    plt.step(np.linspace(-global_args.offset, global_args.offset, num=global_args.offset * 2),
             metagene[1], linewidth=0.5, label=bamname)
    plt.xticks(np.linspace(-global_args.offset, global_args.offset, num=global_args.offset * 2),
               labels, size="xx-small")
    plt.legend()
    plt.savefig(global_args.output_dir + "coverage_term/%s%s.pdf" % (bamname, norm))
    plt.close()
