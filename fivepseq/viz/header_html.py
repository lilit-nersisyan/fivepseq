import os
from fivepseq.logic.structures.fivepseq_counts import FivePSeqCounts

import fivepseq
from fivepseq.logic.structures.annotation import Annotation


def get_version():
    # version_path = os.path.abspath(os.path.join(os.curdir, "version.txt"))
    # if not os.path.exists(version_path):
    #    version_path = os.path.abspath(os.path.join(os.curdir, "fivepseq", "version.txt"))
    # f = open(version_path, 'r')
    # version = f.read()

    return fivepseq.__version__


def write_fivepseq_header(viz_pipeline):
    html_file_path = os.path.join(viz_pipeline.args.o, viz_pipeline.title + ".html")
    html_file = open(html_file_path, "w")

    html_file.write(get_html_body(viz_pipeline))

    html_file.close()


def get_div_logo(w=300, h=100):
    div = """<div>
    <a href="http://pelechanolab.com/software/fivepseq">
    <svg width="%d" height="%d" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 389.27 130.51"><defs><style>.cls-1,.cls-3{fill:#e96656;}.cls-1{fill-rule:evenodd;}.cls-2,.cls-3{font-size:80px;font-family:BookmanOldStyle, Bookman Old Style;}.cls-2{fill:#404040;}</style></defs><title>fivepseq</title><g id="Layer_2" data-name="Layer 2"><g id="Layer_1-2" data-name="Layer 1"><path class="cls-1" d="M51.78,89.51l32-13.81A31,31,0,0,0,58,61.46c-17.42,0-31.54,15.09-31.54,33.69s14.12,33.7,31.54,33.7c12.55,0,23.38-7.83,28.46-19.17Z"/><text class="cls-2" transform="translate(58.5 61.76) rotate(-60.46)">f</text><text class="cls-2" transform="translate(89.53 96.21) rotate(-21.2)">i</text><text class="cls-2" transform="translate(115.15 105.08)">v</text><text class="cls-2" transform="translate(159.03 105.08)">e</text><text class="cls-2" transform="translate(202.85 105.08)">p</text><text class="cls-3" transform="translate(254.68 105.08)">s</text><text class="cls-3" transform="translate(298.49 105.08)">e</text><text class="cls-3" transform="translate(342.31 105.08)">q</text></g></g></svg> 
    </a>
    </div>
    &nbsp;
    &nbsp;""" % (w, h)
    return div


def get_div_title(title):
    div = """
    <h2 style="color:#404040;" align="center">%s</h2>
    """ % title
    return div


def get_div_footer():
    div = """<hr><div>
    These plots are generated with <it>fivepseq</it> version %s. 
    For more information on usage and citation, please visit <a href = "http://pelechanolab.com/software/fivepseq">
    our homepage </a>.
    </div>""" % get_version()
    return div


def get_libsize_table(viz_pipeline):
    table = "<table>" \
            "<tr>" \
            "<th>Sample</th>" \
            "<th>Library size (M)</th>" \
            "</tr>"

    for s in viz_pipeline.samples:
        table += "<tr><td>%s</td><td>%.5f</td></tr>" \
                 % (s, (float(viz_pipeline.lib_size_dict[s]) / 1000000))
    table += "</table>"

    return table


def get_html_body(viz_pipeline):
    args = viz_pipeline.args
    curr_dir = os.path.dirname(__file__)
    with open(os.path.join(curr_dir, "template.html"), "r") as f:
        template = f.read()

        bam_list = ""
        count = 0
        for s in viz_pipeline.samples:
            count += 1
            if count > viz_pipeline.MAX_NUM_SAMPLES:
                bam_list += "<br>" + s + " <font color=\"red\"><i>(not plotted)</i></font>"
            else:
                bam_list += "<br>" + s

        if hasattr(args, "gs") and args.gs is not None:
            input_gene_set = os.path.abspath(args.gs)
            input_geneset_sample_per_gs_report = os.path.join("genesets",
                                                              viz_pipeline.title + "_samples_per_geneset.html")
            input_geneset_genesets_per_sample_report = os.path.join("genesets",
                                                                    viz_pipeline.title + "_genesets_per_sample.html")
        else:
            input_gene_set = "None"
            input_geneset_sample_per_gs_report = "None"
            input_geneset_genesets_per_sample_report = "None"

        if hasattr(args, "gf") and args.gf is not None:
            input_gene_filter = os.path.abspath(args.gf)
        else:
            input_gene_filter = "None"

        if args.no_mask:
            mask = "OFF"
        else:
            if hasattr(args, "codon_mask_size"):
                mask = "by %d positions" % args.codon_mask_size
            else:
                mask = "by %d nts (default)" % FivePSeqCounts.MASK_DIST

        if hasattr(args, "dipeptide_pos"):
            dipeptide_pos = args.dipeptide_pos
        else:
            dipeptide_pos = "%d (default)" % FivePSeqCounts.DIPEPTIDE_POS

        if hasattr(args, "tripeptide_pos"):
            tripeptide_pos = args.tripeptide_pos
        else:
            tripeptide_pos = "%d (default)" % FivePSeqCounts.TRIPEPTIDE_POS

        if viz_pipeline.combine:
            combined_report = os.path.join("main", viz_pipeline.title + "_combined.html")
        else:
            combined_report = "None"

        if hasattr(args, "transcript_type"):
            transcript_type = args.transcript_type
        else:
            transcript_type = Annotation.transcript_type

        result = template.format(svg=get_div_logo(),
                                 version=get_version(),
                                 name=viz_pipeline.title,
                                 title=get_div_title(viz_pipeline.title),
                                 footer=get_div_footer(),
                                 bam_dir=os.path.abspath(args.b),
                                 bam=bam_list,
                                 fa=os.path.abspath(args.g),
                                 gff=os.path.abspath(args.a),
                                 out_dir=args.o,
                                 conflicts=args.conflicts,
                                 gene_set=input_gene_set,
                                 gene_filter=input_gene_filter,
                                 span_size=str(args.span),
                                 transcript_type=transcript_type,
                                 op=args.op,
                                 ds=args.ds,
                                 mask=mask,
                                 tripeptide_pos=tripeptide_pos,
                                 dipeptide_pos=dipeptide_pos,
                                 lib_size_table=get_libsize_table(viz_pipeline),
                                 main_report=os.path.join("main", viz_pipeline.title + "_main.html"),
                                 combined_report=combined_report,
                                 suppl_amino_acid_line_report=os.path.join(
                                     "supplement",
                                     viz_pipeline.title + "_amino_acid_linecharts.html"),
                                 suppl_codon_line_report=os.path.join(
                                     "supplement",
                                     viz_pipeline.title + "_codon_linecharts.html"),
                                 suppl_codon_report=os.path.join(
                                     "supplement",
                                     viz_pipeline.title + "_codon_heatmaps.html"),
                                 suppl_dipeptide_linecharts=os.path.join(
                                     "supplement",
                                     viz_pipeline.title + "_dipeptide_linecharts.html"),
                                 suppl_tripeptide_linecharts=os.path.join(
                                     "supplement",
                                     viz_pipeline.title + "_tripeptide_linecharts.html"),
                                 suppl_dicodon_linecharts=os.path.join(
                                     "supplement",
                                     viz_pipeline.title + "_dicodon_linecharts.html"),
                                 suppl_tricodon_linecharts=os.path.join(
                                     "supplement",
                                     viz_pipeline.title + "_tricodon_linecharts.html"),
                                 differential_heatmaps=os.path.join(
                                     "comparison",
                                     viz_pipeline.title + "_differential_heatmaps.html"
                                 ),
                                 geneset_sample_per_gs_report=input_geneset_sample_per_gs_report,
                                 geneset_genesets_per_sample_report=input_geneset_genesets_per_sample_report)

        return (result)

# if __name__ == '__main__':
#     main()
