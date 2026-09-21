# Notas de setup do ambiente para o `biodenoising`

Resumo dos problemas encontrados ao rodar o comando `biodenoise` neste projeto e das
correções aplicadas no `.venv`. Use este arquivo como referência caso o ambiente
precise ser recriado ou os mesmos erros reapareçam (ex.: após `pip install -U ...`
ou recriação da venv).

## Comando usado

```bash
biodenoise --method biodenoising16k_dns48 \
  --noisy_dir "gravacoes em cativeiro/leem mudancas climaticas/sala 1 09.08.2024" \
  --out_dir "gravacoes em cativeiro/leem mudancas climaticas/sala 1 09.08.2024_biodenoising" \
  --device cuda
```

> Importante: rodar com a venv ativada (`source .venv/bin/activate`) ou via
> `.venv/bin/biodenoise`. Sem ativar, o comando `biodenoise` não é encontrado.

## Problema 1: `ModuleNotFoundError: No module named 'pkg_resources'`

**Causa:** `torchmetrics` (dependência do `biodenoising`) importa `pkg_resources`,
que foi removido do `setuptools` a partir da versão 81.

**Correção:**

```bash
.venv/bin/pip install "setuptools<81"
```

## Problema 2: `module 'torchaudio' has no attribute 'info'`

**Causa:** duas questões combinadas:

1. `torch` estava em uma versão (2.14.0) mais nova do que qualquer `torchaudio`
   disponível no PyPI (máximo 2.11.0) — versões desencontradas.
2. A partir do `torchaudio` 2.9+, a API antiga de I/O (`torchaudio.info`,
   `torchaudio.load`/`save` clássicos) está sendo descontinuada em favor do
   `torchcodec`, e o `biodenoising` (em `denoiser/audio.py`) ainda depende de
   `torchaudio.info(path)`.

**Correção:** fixar um par de versões que ainda expõe `torchaudio.info` e mantém
suporte a CUDA:

```bash
.venv/bin/pip install "torch==2.8.0" "torchaudio==2.8.0"
```

**Verificação:**

```bash
.venv/bin/python -c "
import torch, torchaudio
print(torch.__version__, torchaudio.__version__)
print('has info:', hasattr(torchaudio, 'info'))
print('cuda available:', torch.cuda.is_available())
"
```

Saída esperada: `2.8.0+cu128 2.8.0+cu128`, `has info: True`, `cuda available: True`.

## Checklist rápido para reaplicar tudo do zero

```bash
cd /home/aldemir/Workspace/repositories/chelonian-project
source .venv/bin/activate
pip install "setuptools<81"
pip install "torch==2.8.0" "torchaudio==2.8.0"
```

Depois disso, o comando `biodenoise` deve rodar normalmente (podem aparecer
`UserWarning`s de depreciação do `torchaudio` sobre migração para `torchcodec`,
que são apenas avisos e não impedem a execução).
