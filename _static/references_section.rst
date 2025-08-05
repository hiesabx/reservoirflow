References 📦
--------------
.. usage in rst files:
    .. include:: /_static/comments_section.rst

.. usage in ipynb files: 
    1. you need to remove .. raw:: directive
    ```{include} /_static/comments_section.rst
    :heading-offset: 1
    ```
    1. with .. raw:: directive but you need to add a header.
    ```{eval-rst}
    .. include:: /_static/comments_section.rst
        :start-line: 3
    ```

.. footbibliography::
