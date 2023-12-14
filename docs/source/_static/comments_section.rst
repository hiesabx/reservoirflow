Comments 💬
-----------
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


Feel free to make a comment, ask a question, or share your opinion about this specific content. 
Please keep in mind the `Commenting Guidelines ⚖ </community/commenting_guidelines.html>`_.

.. raw:: html
    :class: only-dark

    <script 
        type="text/javascript"
        src="https://utteranc.es/client.js"
        async="async"
        repo="zakgrin/reservoirflow_utterances"
        issue-term="pathname"
        theme="github-dark"
        label="comments 💬"
        crossorigin="anonymous"
    />
    </script>

.. raw:: html
    :class: only-light

    <script 
        type="text/javascript"
        src="https://utteranc.es/client.js"
        async="async"
        repo="zakgrin/reservoirflow_utterances"
        issue-term="pathname"
        theme="github-light"
        label="comments 💬"
        crossorigin="anonymous"
    />
    </script>
