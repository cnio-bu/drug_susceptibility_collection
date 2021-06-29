rule annotate_cell_lines:
    input:
        sample_info=datasets.loc['ccle_sample_info', 'directory'],
        celligner_lines_data=datasets.loc['celligner_data', 'directory']
    output:
        cell_lines_annotation=f'{results}/annotation/cell_line_annotation.csv'
    threads: 1
    resources:
        mem_mb=4096
    conda:
        '../envs/file_manipulation.yaml'
    script:
        '../scripts/common_annotate_lines.py'